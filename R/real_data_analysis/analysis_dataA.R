# ----------load data 
scNGP.C1 <- readRDS(choose.files()) # choose data A
scNGP.C1 <- scNGP.C1$filtered  # the data contains preprocessed SingleCellExperiment objects

# filter genes
treatment <- scNGP.C1$treatment
table(treatment)
keep.gene <- apply(counts(scNGP.C1), 1, function(y){
  #all(tapply(y, treatment, function(x) sum(x>0)>5 & sum(y)>=5))
  (sum(y>0)>=5)  &  sum(y)>=5
})
table(keep.gene)

scNGP.C1_filtered <- scNGP.C1[keep.gene, ]
scNGP.C1_filtered 
summary(rowSums(counts(scNGP.C1_filtered)))
summary(rowSums(counts(scNGP.C1_filtered)>0))

# extra cell level info
scNGP.C1_filtered$treatment <- ifelse(scNGP.C1_filtered$treatment == "nutlin", 2, 1)
scNGP.C1_filtered$lLS       <- log(colSums(counts(scNGP.C1_filtered)))


#  ---------- run DGE analysis
library(SingleCellExperiment)

source(".../tools_wrap_functions.R")
 
res.pim <- PIMSeq(SCExp =  scNGP.C1_filtered, condition =  "treatment",
                        nuisance.vars =  "lLS", assay.name = "counts",
                  BPPARAM = BiocParallel::SnowParam())  # run PIMSeq with log-library size in the model

res.SAMSeq <- run_SAMseq(SCEdata = scNGP.C1_filtered, condition.name = "treatment",
                       covariates = NULL, expression.unit = "count", prop.zero = 1) # SAMseq 
res.MAST <- run_MAST(SCEdata = scNGP.C1_filtered, condition.name = "treatment",
                         covariates = NULL, expression.unit = "CPM", prop.zero = 1) # MAST with cellular detection rate in the model
res.edgeR_Zinger <- run_edgeR_Zinger(SCEdata = scNGP.C1_filtered, condition.name = "treatment",
                     covariates = NULL, expression.unit = "count", prop.zero = 1)  #edgeR +Zinger
res.DESeq2_Zinger <- run_DESeq2_Zinger(SCEdata = scNGP.C1_filtered, condition.name = "treatment",
                                     covariates = NULL, expression.unit = "count", prop.zero = 1) #edgeR + DESeq2

saveRDS(res.pim, ".../scNGP_C1_analysis/res.pim.rds")
saveRDS(res.SAMSeq, ".../scNGP_C1_analysis/res.SAMSeq.rds")
saveRDS(res.MAST, ".../scNGP_C1_analysis/res.MAST.rds")
saveRDS(res.edgeR_Zinger, ".../scNGP_C1_analysis/res.edgeR_Zinger.rds")
saveRDS(res.DESeq2_Zinger, ".../scNGP_C1_analysis/res.DESeq2_Zinger.rds")


#  ----------  UPset plot (Figure 4 -- panel A) -------------------------------------
library(UpSetR)
DE_class.pim    <- as.numeric(res.pim$test.contrasts$p.adjusted < 0.05 )
DE_class.SAMSeq <- as.numeric(res.SAMSeq$df$qval < 0.05 & res.SAMSeq$df$FoldChange>0)
DE_class.MAST   <- as.numeric(res.MAST$df$qval < 0.05)
DE_class.edgeR  <- as.numeric(res.edgeR_Zinger$df$qval < 0.05)
DE_class.DESeq2 <- as.numeric(res.DESeq2_Zinger$df$padj < 0.05)

 
DE.class <- as.data.frame(cbind(DE_class.pim, DE_class.SAMSeq,  
                                DE_class.MAST, DE_class.edgeR, DE_class.DESeq2))
rownames(DE.class) <- rownames(scNGP.C1_filtered) 
colnames(DE.class) <- c("PIM",  "SAMSeq", "MAST", "edgeR+Zinger", "DESeq2+Zinger")
DE.class <- DE.class[complete.cases(DE.class),]
head(DE.class)
colSums(DE.class)

saveRDS(DE.class, ".../scNGP_C1_analysis/DE.class.rds")

png(filename = ".../scNGP_C1_analysis/upSetPlot.png", width=15, 
    height=10, res = 800, units = "cm")
upset(as.data.frame(DE.class), order.by = c("freq"), decreasing = c(TRUE),
              mainbar.y.label="", sets.x.label="", att.pos="top", cutoff=1)
dev.off()

# ---------- Unque set of DE genes   (Figure S7 --supllementary file)  ------------------------
SAMSeq.only.DE <- rownames(DE.class[DE.class$SAMSeq==1 & DE.class$PIM==0 & DE.class$MAST==0 & 
                             DE.class$`edgeR+Zinger`==0 & DE.class$`DESeq2+Zinger`==0, ])
PIM.only.DE    <- rownames(DE.class[DE.class$SAMSeq==0 & DE.class$PIM==1 & DE.class$MAST==0 & 
                             DE.class$`edgeR+Zinger`==0 & DE.class$`DESeq2+Zinger`==0, ])
MAST.only.DE   <- rownames(DE.class[DE.class$SAMSeq==0 & DE.class$PIM==0 & DE.class$MAST==1 & 
                             DE.class$`edgeR+Zinger`==0 & DE.class$`DESeq2+Zinger`==0, ])
edgeR.only.DE  <- rownames(DE.class[DE.class$SAMSeq==0 & DE.class$PIM==0 & DE.class$MAST==0 & 
                          DE.class$`edgeR+Zinger`==1 & DE.class$`DESeq2+Zinger`==0, ])
DESeq2.only.DE <- rownames(DE.class[DE.class$SAMSeq==0 & DE.class$PIM==0 & DE.class$MAST==0 & 
                            DE.class$`edgeR+Zinger`==0 & DE.class$`DESeq2+Zinger`==1, ])
unique.DE.genes <- list(PIM=PIM.only.DE, SMASeq=SAMSeq.only.DE, MAST=MAST.only.DE,
                        edgeR_Zinger=edgeR.only.DE, DESeq2_Zinger=DESeq2.only.DE)

cpm.mat  <- edgeR::cpm(counts(scNGP.C1_filtered))
var.CPM  <- rowVars(cpm.mat)
foldChange.CPM  <- apply(cpm.mat,1,function(x){
  abs(diff(tapply(log2(x+1), scNGP.C1_filtered$treatment, mean)))
})
extrem.var.genes             <- rownames(cpm.mat)[var.CPM>10000]
small.foldChange.genes       <- rownames(cpm.mat)[foldChange.CPM<0.5]

extreme.DE.genes_var <- sapply(unique.DE.genes, function(g){
  sum(g %in% extrem.var.genes, na.rm = TRUE)
})

small.FC.DE.genes <- sapply(unique.DE.genes, function(g){
  sum(g %in% small.foldChange.genes, na.rm = TRUE)
})


extreme.DE.genes <- data.frame(extreme_var=extreme.DE.genes_var, 
                               small_FC=small.FC.DE.genes)
extreme.DE.genes$tool <- rownames(extreme.DE.genes)
extreme.DE.genes <- reshape2::melt(extreme.DE.genes, id.vars="tool")
extreme.DE.genes$Characteristics <- as.character(extreme.DE.genes$variable)
extreme.DE.genes[extreme.DE.genes$variable=="extreme_var", "Characteristics"]<-"extreme variability"
extreme.DE.genes[extreme.DE.genes$variable=="small_FC", "Characteristics"]   <- "logFC<0.5"

library(ggplot2)
png(filename = ".../scNGP_C1_analysis/DE_genes_xstics.png", width=15, 
    height=10, res = 800, units = "cm")
ggplot(extreme.DE.genes, aes(x=tool, y=value+1, group=Characteristics, fill=Characteristics))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x=NULL, y="number of unique DE genes (at 5% FDR)")+
  coord_flip()+
  theme_bw()+
  theme(legend.position = c(0.7, 0.3),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17))
dev.off()


# ------------ distribution of difference in zero fraction among DE genes (Figure 4 panle B) -----------------
frac.zero.diff <- apply(counts(scNGP.C1_filtered), 1, function(y){
  abs(diff(tapply(y, scNGP.C1_filtered$treatment, function(x) mean(x==0))))
})
head(frac.zero.diff)
frac.zero.diff.pim     <- frac.zero.diff[res.pim$test.contrasts$p.adjusted<.05]
frac.zero.diff.SAMSeq  <- frac.zero.diff[res.SAMSeq$df$qval<.05]
frac.zero.diff.MAST    <- frac.zero.diff[res.MAST$df$qval<.05]
frac.zero.diff.edgeR   <- frac.zero.diff[res.edgeR_Zinger$df$qval<.05]
frac.zero.diff.DESeq2  <- frac.zero.diff[res.DESeq2_Zinger$df$padj<.05]

par(mfrow=c(1,5))
boxplot(frac.zero.diff.pim, main="PIM")
boxplot(frac.zero.diff.SAMSeq, main="SAMSeq")
boxplot(frac.zero.diff.MAST, main="MAST")
boxplot(frac.zero.diff.edgeR, main="edgeR+Zinger")
boxplot(frac.zero.diff.DESeq2, main="DESeq2+Zinger")
par(mfrow=c(1,1))



# -------------   Explore unique set of DE genes by each tool: Beeswarm plot (Figure S6 --supllementary file)  ------------------------
C1.data.cpm <- log2(edgeR::cpm(counts(scNGP.C1_filtered))+1)
library(beeswarm)
set.seed(234)
grp <- scNGP.C1_filtered$treatment  

png(".../scNGP_C1_analysis/beeswarmPlot.png", width = 30, height = 25,
    units = "cm", res = 500)
par(mfrow=c(5, 5), mai=c(0.3, 0.3, 0.5, 0.1))
genes_PIM_only <- rownames(DE.class[DE.class$PIM==1 & DE.class$SAMSeq==0 & DE.class$MAST==0 & DE.class$`edgeR+Zinger`==0 & DE.class$`DESeq2+Zinger`==0,])
genes_PIM_only <- sample(genes_PIM_only, 5) 
for(gg in genes_PIM_only){
  count_only <- C1.data.cpm[gg, ] 
  beeswarm(count_only ~ grp, 
           method ="center", corral="gutter", bg = "#00000050",
           pch = 16,  xlab = "", ylab = "log2(CPM+1)", 
           labels=c("vehicle", "nutlin"), col=c(1,2),
           main=paste0(gg, "\n", "FDR =", 
                       round(res.pim$df$p.adjusted[res.pim$df$ID==gg], 4),
                       "   PI=", round(res.pim$df$PI[res.pim$df$ID==gg], 2)))
} 

genes_SAMSeq_only <- rownames(DE.class[DE.class$PIM==0 & DE.class$SAMSeq==1 & DE.class$MAST==0 & DE.class$`edgeR+Zinger`==0 & DE.class$`DESeq2+Zinger`==0,])
genes_SAMSeq_only <- sample(genes_SAMSeq_only, 5)
for(gg in genes_SAMSeq_only){
  count_only <- C1.data.cpm[gg, ] 
  beeswarm(count_only ~ grp, 
           method ="center", corral="gutter", bg = "#00000050",
           pch = 16,  xlab = "", ylab = "log2(CPM+1)", 
           labels=c("vehicle", "nutlin"),  col=c(1,2),
           main=paste0(gg, "\n", "FDR =", 
                       round(res.SAMSeq$df[gg, "qval"], 4)))
}

genes_MAST_only <- rownames(DE.class[DE.class$PIM==0 & DE.class$SAMSeq==0 & DE.class$MAST==1 & DE.class$`edgeR+Zinger`==0 & DE.class$`DESeq2+Zinger`==0,])
genes_MAST_only <- sample(genes_MAST_only, 5)
for(gg in genes_MAST_only){
  count_only <- C1.data.cpm[gg, ] 
  beeswarm(count_only ~ grp, 
           method ="center", corral="gutter", bg = "#00000050",
           pch = 16,  xlab = "", ylab = "log2(CPM+1)", 
           labels=c("vehicle", "nutlin"), col=c(1,2),
           main=paste0(gg, "\n", "FDR =", 
                       round(res.MAST$df[gg, "qval"], 4)))
}


genes_edgeR_DESeq2_only <- rownames(DE.class[DE.class$PIM==0 & DE.class$SAMSeq==0 & DE.class$MAST==0 & (DE.class$`edgeR+Zinger`==1 | DE.class$`DESeq2+Zinger`==1),])
genes_edgeR_DESeq2_only <- sample(genes_edgeR_DESeq2_only, 5)
for(gg in genes_edgeR_DESeq2_only){
  count_only <- C1.data.cpm[gg, ] 
  beeswarm(count_only ~ grp, 
           method ="center", corral="gutter", bg = "#00000050",
           pch = 16,  xlab = "", ylab = "log2(CPM+1)", 
           labels=c("vehicle", "nutlin"), col=c(1,2),
           main=paste0(gg, "\n", "FDR =", 
                       round(res.edgeR_Zinger$df[gg, "qval"], 4)))
}


common_genes <- rownames(DE.class[DE.class$PIM==1 & DE.class$SAMSeq==1 & DE.class$MAST==1 & DE.class$`edgeR+Zinger`==1 & DE.class$`DESeq2+Zinger`==1,])
common_genes <- sample(common_genes, 5)
#par(mfrow=c(1, 5))
for(gg in common_genes){
  count_only <- C1.data.cpm[gg, ] 
  beeswarm(count_only ~ grp, 
           method ="center", corral="gutter", bg = "#00000050",
           pch = 16,  xlab = "", ylab = "log2(CPM+1)", 
           labels=c("vehicle", "nutlin"), col=c(1,2),
           main=paste0(gg, "\n", "(PIM) FDR =", 
                       round(res.pim$df$p.adjusted[res.pim$df$ID==gg], 4),
                       "   PI=", round(res.pim$df$PI[res.pim$df$ID==gg], 2)))
}
dev.off()



#------------------ Fold change estimate comparison --- (Figure 5 panel A, B and C)  --------------------------
plotPIrank2 <- function(pim.res, contrast, add.PI.interval=TRUE, FDR.thrld=0.05,  
         PI.thrld=c(0.4, 0.6), xlab="gene rank", ylab="probabilstic index", 
         add.legend=TRUE, cols.seg=NULL, las=1, cex.lab=1.25, ylim=c(0,1), 
         col=NULL, pch=20, cex=0.5, lty.PI=3, lwd.PI=1.25, col.PI="blue",
         legend.pos= "bottomright", pt.cex.legend = 1.25, 
         cex.legend=1, ...){
  
  expit  <- function(x){exp(x)/(1+exp(x))}
  
  df <- testPIMcontrast(pim.res, contrasts = contrast)
  
  #df <- merge(pim.contrast, pim.res$all.coefficients, by="ID") 
  
  if(add.PI.interval){
    df$lPI <- expit(df$contrast-1.96*df$std.error)
    df$uPI <- expit(df$contrast+1.96*df$std.error)
  }
  df  <- df[order(df$PI), ] 
  
  if(is.null(cols.seg)) cols.seg <- ifelse(df$p.value>=FDR.thrld, "gray", "rosybrown2")
  
  plot(1:nrow(df), df$PI, type="n", xlab=xlab,
       ylab = ylab, las=las, cex.lab=cex.lab, ylim=ylim, ...)
  if(add.PI.interval){
    segments(1:nrow(df), df$lPI, 1:nrow(df), df$uPI, col=cols.seg, ...) 
  }
  
  if(is.null(col)){
    col <- factor(df$p.adjusted<FDR.thrld)
  }
  points(1:nrow(df), df$PI, col=col, pch=pch, cex=cex, ...)
  
  abline(h=PI.thrld, col=col.PI, lty=lty.PI, lwd=lwd.PI) 
  
  if(add.legend){
    legend("bottomright", c(paste("FDR>=", FDR.thrld), paste("FDR<", FDR.thrld)), 
           col = unique(col),  pch = pch, cex=cex.legend, pt.cex = pt.cex.legend)
  } 
}
library(RColorBrewer)

png(filename = ".../scNGP_C1_analysis/PIranking1.png", width=16, height=10, res = 300, units = "cm") 
par(oma=c(0.1,0.1,0.1,0.1), mar=c(4, 4, 0.1, 0.1))
layout(matrix(c(1,1,2,3), 2,2, byrow = FALSE), width=c(2, 1))
plotPIrank2(res.pim$res, contrast = c(1), axes = FALSE, ann = TRUE, PI.thrld =  0.5)
axis(1, at= c(1, 4000, 8000, 12000, 16525), c(1, 4000, 8000, 12000, 16525), col="gray", cex.axis=0.75)
axis(2, at=seq(0, 1,  0.25), seq(0, 1, 0.25), col="gray", las=1, cex.axis=0.75)
text(4500, 0.1, paste("#down-reg. genes =", sum(res.pim$df$p.adjusted<0.05 & res.pim$df$PI<0.5)),
     col="blue")
text(12000, 0.9, paste("#up-reg. genes =", sum(res.pim$df$p.adjusted<0.05 & res.pim$df$PI>0.5)),
     col="blue")
box("figure", col="gray") 


res.PIM_DESeq2 <- merge(res.pim$df, as.data.frame(res.DESeq2_Zinger$df), by.x="ID", by.y="row.names")
par(mar=c(3, 3, 0.1, 0.1),bty="n", fg="gray", cex=0.75)
#Lab.palette <- colorRampPalette(c("white", "red"), space = "Lab")
smoothScatter(res.PIM_DESeq2$PI, -1*res.PIM_DESeq2$log2FoldChange, xlim=c(0, 1), 
              nrpoints =150, las=1, cex.axis=0.75)
title(ylab="log2 FC", xlab="probabilistc index", line=2)


par(mar=c(3, 0.1, 0.1, 0.1))
hist(res.pim$df$PI, nclass = 50, col="gray50", border = "cyan", xlab="", 
     xlim = c(0, 1), las=1, main = "", ylab="", cex.axis=0.75, yaxt='n')
title(ylab="", line=0)
title(xlab="probabilistc index", line=2)
box("figure", col="gray")
dev.off()




# --------------  Gene set enrichment analysis (Figure 4 panel C) ---------------------
library(fgsea)
library(hgu95av2.db) 

gene.set <- read.table(".../Top_116_genes_p53_activated.txt",  sep="\t", header = TRUE)
head(gene.set)
ensembl <- mapIds(hgu95av2.db, keys=as.character(gene.set$Gene.symbol),
                  keytype="SYMBOL",  column="ENSEMBL")
gene.set$ENSEMBL <- ensembl
head(gene.set) 
dim(gene.set)   
pathways <- as.character(gene.set$ENSEMBL) 
pathways <- pathways[!is.na(pathways)]

res.SAMSeq$df$`Score(d)`[is.na(res.SAMSeq$df$`Score(d)`)] <- 0

stat.set.PIM      <- res.pim$res$all.coefficients$`beta:treatment`/res.pim$res$all.coefficients$`SE:treatment`
stat.set.SAMSeq   <- res.SAMSeq$df$`Score(d)` 
stat.set.MAST     <- res.MAST$df$test.stat #C1.results.list2$res.MAST_A$pval
stat.set.edgeR    <- qnorm(1- res.edgeR_Zinger$tt$table$PValue/2)*sign(res.edgeR_Zinger$tt$table$logFC)
res.DESeq2_Zinger$df$pvalue[is.na(res.DESeq2_Zinger$df$pvalue)] <- 1
stat.set.DESeq2  <- qnorm(1-res.DESeq2_Zinger$df$pvalue/2)*sign(-1*res.DESeq2_Zinger$df$log2FoldChange)


names(stat.set.PIM)    <- res.pim$res$all.coefficients$ID
names(stat.set.SAMSeq) <- rownames(res.SAMSeq$df) #C1.results.list2$res.SAMSeq_A$ID
names(stat.set.MAST)   <- rownames(res.MAST$df) #C1.results.list2$res.MAST_A$ID
names(stat.set.edgeR)  <- rownames(res.edgeR_Zinger$tt$table)
names(stat.set.DESeq2) <- rownames(res.DESeq2_Zinger$df)
 

fgseaRes.PIM   <- fgsea(list(p58=pathways), stat.set.PIM, minSize=15, 
                        maxSize=500, nperm=10000)
fgseaRes.SAMSeq   <- fgsea(list(p58=pathways), stat.set.SAMSeq, minSize=15, 
                           maxSize=500, nperm=10000)
fgseaRes.MAST   <- fgsea(list(p58=pathways), stat.set.MAST, minSize=15, 
                         maxSize=500, nperm=10000)
fgseaRes.edgeR   <- fgsea(list(p58=pathways), stat.set.edgeR, minSize=15, 
                          maxSize=500, nperm=10000)
fgseaRes.DESeq2 <- fgsea(list(p58=pathways), stat.set.DESeq2, 
                         minSize=15, maxSize=500, nperm=10000)

score.df <- rbind(fgseaRes.PIM, fgseaRes.SAMSeq, fgseaRes.MAST, fgseaRes.edgeR, fgseaRes.DESeq2)
score.df$tools <- c("PIM", "SAMSeq", "MAST", "edgeR+Zinger", "DESeq2+Zinger")
head(score.df) 

score.df2 <- score.df[, c(2,4,5,7, 9)]
score.df3 <- reshape2::melt(score.df2[, c(2, 3, 5)])
colnames(score.df3) <- c("tools", "type", "score")
ggplot(score.df3, aes(x=tools, y=score, fill=type))+
  geom_bar(stat = "identity", position=position_dodge(), width=0.25)+
  labs(x=NULL)+
  coord_flip()+
  theme_classic()

score.df3 <- score.df2[order(score.df2$NES),]

png(filename = ".../scNGP_C1_analysis/GSEAsummary.png", 
    width = 10, height = 5, units = "cm", res = 300)
par(mar=c(2, 7, 0.1, 0.1 ))
plot(score.df3$NES, 1:nrow(score.df3),  xlim = c(0.5, 4), ylim=c(1, 5.25),
     axes = FALSE, ylab="", xlab="", pch=19, cex=1.5)
segments(rep(0, nrow(score.df3)), 1:nrow(score.df3),  score.df3$NES, 1:nrow(score.df3), 
         col="gray", lwd=1.25)
points(score.df3$ES, (1:nrow(score.df3))+0.1,  col="dodgerblue3", pch=18, cex=1.5)
segments(rep(0, nrow(score.df3)), (1:nrow(score.df3))+0.1, score.df3$ES,  (1:nrow(score.df3))+0.1, 
         col="cadetblue2", lwd=1.25)
axis(2, 1:5, score.df3$tools, las=1)
axis(1, seq(0.5, 4, 1), las=1)

legend("bottomright",c("NES", "ES"), col=c("black", "dodgerblue3"), pch = c(19, 18), lty=1,
       cex = 0.75)
box("figure", col="gray")
dev.off()
abline(v=0.5, h=0, lty=3, col="gray60")
box("figure", col="gray")


par(mar=c(3, 0.1, 0.1, 0.1))
hist(res.pim$df$PI, nclass = 50, col="gray50", border = "cyan", xlab="", 
     xlim = c(0, 1), las=1, main = "", ylab="", cex.axis=0.75, yaxt='n')
title(ylab="", line=0)
title(xlab="probabilistc index", line=2)
box("figure", col="gray")
dev.off()

#####################################################################################################################
#####################################################################################################################           
           
# ------------- Analysis of NGP Chromium data -------------
# load data
scNGP.10x <- readRDS(choose.files())  # choose dataB
scNGP.10x <- scNGP.10x$filtered
summary(rowSums(counts(scNGP.10x)))
summary(rowSums(counts(scNGP.10x)>0))

# filter genes
treatment <- scNGP.10x2$treatment
table(treatment)
keep.gene <- apply(counts(scNGP.10x2), 1, function(y){
  #all(tapply(y, treatment, function(x) sum(x>0)>3))
  (sum(y>0)>=3)  &  sum(y)>=3
})
table(keep.gene)


scNGP.10x2_filtered <- scNGP.10x2[keep.gene, ]
scNGP.10x2_filtered 
summary(rowSums(counts(scNGP.10x2_filtered)))
summary(rowSums(counts(scNGP.10x2_filtered)>0))


# ---------------run DGE  analysis --------------------------------------

source(".../tools_wrap_functions.R")
scNGP.10x2_filtered$lLS       <- log(colSums(counts(scNGP.10x2_filtered)))       
res.pim <- run_PIMseq(SCEdata = scNGP.10x2_filtered, condition.name = "treatment",
                      covariates = "lLS", prop.zero = 1, coxph.aprox = TRUE)
res.SAMSeq <- run_SAMseq(SCEdata = scNGP.10x2_filtered, condition.name = "treatment",
                         covariates = NULL, expression.unit = "count", prop.zero = 1)  # SAMSeq is failed
res.MAST <- run_MAST(SCEdata = scNGP.10x2_filtered, condition.name = "treatment",
                     covariates = NULL, expression.unit = "CPM", prop.zero = 1)
res.edgeR_Zinger <- run_edgeR_Zinger(SCEdata = scNGP.10x2_filtered, condition.name = "treatment",
                                     covariates = NULL, expression.unit = "count", prop.zero = 1)
res.DESeq2_Zinger <- run_DESeq2_Zinger(SCEdata = scNGP.10x2_filtered, condition.name = "treatment",
                                       covariates = NULL, expression.unit = "count", prop.zero = 1)

saveRDS(res.pim, ".../scNGP_10x_analysis/res.pim.rds")
#saveRDS(res.SAMSeq, ".../scNGP_10x_analysis/res.SAMSeq.rds")
saveRDS(res.MAST, ".../scNGP_10x_analysis/res.MAST.rds")
saveRDS(res.edgeR_Zinger, ".../scNGP_10x_analysis/res.edgeR_Zinger.rds")
saveRDS(res.DESeq2_Zinger, "..../scNGP_10x_analysis/res.DESeq2_Zinger.rds")

           
           
# ---------------  cross data agreement (Figure 4 -- panel D)----------------
DE_class.pim    <- as.numeric(res.pim$df$p.adjusted < 0.05 & (res.pim$df$PI<=0.4 | res.pim$df$PI>=0.6))
DE_class.MAST   <- as.numeric(res.MAST$df$qval < 0.05 & abs(res.edgeR_Zinger$tt$table$logFC)>0.6)
DE_class.edgeR  <- as.numeric(res.edgeR_Zinger$df$qval < 0.05 & abs(res.edgeR_Zinger$tt$table$logFC)>0.6)
DE_class.DESeq2 <- as.numeric(res.DESeq2_Zinger$df$padj < 0.05 & abs(res.DESeq2_Zinger$df$log2FoldChange)>0.6)

DE.class_10x <- as.data.frame(cbind(DE_class.pim,  DE_class.MAST, 
                                    DE_class.edgeR, DE_class.DESeq2))
rownames(DE.class_10x) <- rownames(scNGP.10x_filtered) 
colnames(DE.class_10x) <- c("PIM",   "MAST", "edgeR+Zinger", "DESeq2+Zinger")
DE.class_10xs <- DE.class_10x[complete.cases(DE.class_10x),]
head(DE.class_10x)
colSums(DE.class_10x)

DE.class_C1 <- readRDS(".../scNGP_C1_analysis/DE.class.rds")
head(DE.class_C1)
colSums(DE.class_C1)

common.genes <- intersect(rownames(DE.class_C1), rownames(DE.class_10x))
length(common.genes)

DE.class_C12 <- DE.class_C1[common.genes, -2]
DE.class_10x2 <- DE.class_10x[common.genes,]

common.DE.genes <- t(sapply(common.genes, function(g){
  v <- numeric(4)
  for(i in 1:4){
    v[i] <- as.numeric((DE.class_C12[g, i])==1 & (DE.class_10x2[g, i]==1))
  }
  v
}))
head(common.DE.genes)
colSums(common.DE.genes)

cross.data.res <- data.frame(C1=colSums(DE.class_C1[, -2]),
                             Chromium=colSums(DE.class_10x),
                             common=colSums(common.DE.genes),
                             tool=c("PIM", "MAST", "edgeR+Zinger", "DESeq2+Zinger"))
cross.data.res
cross.data.res2 <- reshape2::melt(cross.data.res, id.vars="tool")
library(ggplot2)
png(filename = "manuscript/V3/supplementary file/scNGP_10x_analysis/cross_data.png", width=5, height=7, 
    res = 800, units = "cm")
ggplot(cross.data.res2, aes(x=tool, y=value, group=variable, fill=variable))+
  geom_bar(stat="identity", position=position_dodge(), width = 0.6, size=0.2, colour="black") +
  #geom_point(size=3, aes(x=tool, y=value, fill=variable), position=position_dodge()) +
  labs(y="#DE genes (at 5% FDR)", x=NULL)+
  #scale_fill_brewer(palette="Paired") + 
  scale_fill_manual(values=c("deepskyblue4", '#999999', 'darkorange1'))+
  theme_minimal() + 
  theme(axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8, angle = 45, hjust = 1),
        axis.text.y = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.25, 'cm'),
        legend.spacing.x = unit(0.15, 'cm')) 
dev.off()
