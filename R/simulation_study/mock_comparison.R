# ----------------------- Mock comparison with n=21 biological replicates --------------
library(SingleCellExperiment)
scNGP.polyA <- readRDS("Data/scNGP_polyA.rds")
colData(scNGP.polyA)

scNGP.polyA_vehicle <- scNGP.polyA[, scNGP.polyA$characteristics..treatment=="vehicle"]
scNGP.polyA_vehicle <- scNGP.polyA_vehicle[rowSums(counts(scNGP.polyA_vehicle)>0)>5,]

set.seed(6412)
keep <- sample(nrow(scNGP.polyA_vehicle), 4000)
scNGP.polyA_vehicle <- scNGP.polyA_vehicle[keep, ]  

itr <- 100
X <- rep(0:1, each=ncol(scNGP.polyA_vehicle)/2)
library(permute)
mock.X <- sapply(1:itr, function(i){
  set.seed(264+i)
  X[shuffle(X)]
})

source("tools_wrap_functions.R")
# 
mock.res.pim.list <- lapply(1:itr, function(i){
  print(i)
  mock.group <- mock.X[, i]
  scData            <- scNGP.polyA_vehicle
  scData$mock.group <- mock.group
  scData$lLS        <- log(colSums(counts(scData)))
  
  res <- run_PIMseq(SCEdata = scData, condition.name = "mock.group", covariates = "lLS",
                    prop.zero = 0.99, coxph.aprox = FALSE)
  res$df
})
saveRDS(mock.res.pim.list, ".../mock_comparison/dataC1/mock.res.pim.list.rds")

mock.res.pim.list_CoxPH <- lapply(1:itr, function(i){
  print(i)
  mock.group <- mock.X[, i]
  scData            <- scNGP.polyA_vehicle
  scData$mock.group <- mock.group
  scData$lLS        <- log(colSums(counts(scData)))
  
  res <- run_PIMseq(SCEdata = scData, condition.name = "mock.group", covariates = "lLS",
                    prop.zero = 0.99, coxph.aprox = TRUE)
  res$df
})
saveRDS(mock.res.pim.list_CoxPH, ".../mock_comparison/dataC1/mock.res.pim.list_CoxPH.rds")

mock.res.MAST.list <- lapply(1:itr, function(i){
  print(i)
  mock.group <- mock.X[, i]
  scData            <- scNGP.polyA_vehicle
  scData$mock.group <- mock.group
  
  res <- run_MAST(SCEdata = scData, condition.name = "mock.group", covariates = NULL,
                  expression.unit = "CPM", prop.zero = 0.99)
  res$df
})
saveRDS(mock.res.MAST.list, ".../mock_comparison/dataC1/mock.res.MAST.list.rds")


mock.res.edgeR_Zinger.list <- lapply(1:itr, function(i){
  print(i)
  mock.group <- mock.X[, i]
  scData            <- scNGP.polyA_vehicle
  scData$mock.group <- mock.group
  
  res <- run_edgeR_Zinger(SCEdata = scData, condition.name = "mock.group", covariates = NULL,
                          expression.unit = "count", prop.zero = 0.99)
  res$df
})
saveRDS(mock.res.edgeR_Zinger.list, ".../mock_comparison/dataC1/mock.res.edgeR_Zinger.list.rds")


mock.res.DESeq2_Zinger.list <- lapply(1:itr, function(i){
  print(i)
  mock.group <- mock.X[, i]
  scData            <- scNGP.polyA_vehicle
  scData$Group <- mock.group
  
  res <- run_DESeq2_Zinger(SCEdata = scData, condition.name = "Group", covariates = NULL,
                           expression.unit = "count", prop.zero = 0.99)
  res$df
})
saveRDS(mock.res.DESeq2_Zinger.list, ".../mock_comparison/dataC1/mock.res.DESeq2_Zinger.list.rds")


 
boxplot(list(PIM = rowMeans(sapply(mock.res.pim.list, function(x) as.numeric(x$p.value<0.05))),
             PIM_CoxPH = rowMeans(sapply(mock.res.pim.list2, function(x) as.numeric(x$p.value<0.05))),
             MAST = rowMeans(sapply(mock.res.MAST.list, function(x) as.numeric(x$pval<0.05))),
             DESeq2_Zinger=rowMeans(sapply(mock.res.DESeq2_Zinger.list, function(x) as.numeric(x$pvalue<0.05))),
             edgeR_Zinger=rowMeans(sapply(mock.res.edgeR_Zinger.list, function(x) as.numeric(x$pval<0.05)))),
        ylim=c(0, 1))
abline(h=0.05, col=2, lty=2)









# ----------------------- Mock comparison with n=50 biological replicates --------------
library(SingleCellExperiment)
scNGP.10x <- readRDS("Data/scNGP10x_data.rds")
colData(scNGP.10x)


scNGP.10x_vehicle <- scNGP.10x[, scNGP.10x$treatment=="control"]
scNGP.10x_vehicle <- scNGP.10x_vehicle[rowSums(counts(scNGP.10x_vehicle)>0)>5,]
set.seed(6713)
scNGP.10x_vehicle <- scNGP.10x_vehicle[sample(nrow(scNGP.10x_vehicle), 4000),]
  
  
itr <- 100
X <- rep(0:1, each=50)
library(permute)
mock.X <- sapply(1:itr, function(i){
  set.seed(264+i)
  sample(ncol(scNGP.10x_vehicle), 100)
})


source("tools_wrap_functions.R")
# 
mock.res.pim.list <- lapply(1:itr, function(i){
  print(i) 
  scData            <- scNGP.10x_vehicle[, mock.X[, i]]
  scData$mock.group <- X
  scData$lLS        <- log(colSums(counts(scData)))

  res <- run_PIMseq(SCEdata = scData, condition.name = "mock.group", covariates = "lLS",
             prop.zero = 0.99, coxph.aprox = FALSE,  ncores = 7)
  res$df
})
saveRDS(mock.res.pim.list, ".../mock_comparison/data10x/mock.res.pim.list.rds")

mock.res.pim.list_CoxPH <- lapply(1:itr, function(i){
  print(i) 
  scData            <- scNGP.10x_vehicle[, mock.X[, i]]
  scData$mock.group <- X
  scData$lLS        <- log(colSums(counts(scData)))
  
  res <- run_PIMseq(SCEdata = scData, condition.name = "mock.group", covariates = "lLS",
                    prop.zero = 0.99, coxph.aprox = TRUE)
  res$df
})
saveRDS(mock.res.pim.list_CoxPH, ".../mock_comparison/data10x/mock.res.pim.list_CoxPh.rds")

mock.res.MAST.list <- lapply(1:itr, function(i){
  print(i)
  scData            <- scNGP.10x_vehicle[, mock.X[, i]]
  scData$mock.group <- X

  res <- run_MAST(SCEdata = scData, condition.name = "mock.group", covariates = NULL,
                  expression.unit = "CPM", prop.zero = 0.99)
  res$df
})
saveRDS(mock.res.MAST.list, ".../mock_comparison/data10x/mock.res.MAST.list.rds")
 
mock.res.edgeR_Zinger.list <- lapply(1:itr, function(i){
  print(i)
  scData            <- scNGP.10x_vehicle[, mock.X[, i]]
  scData$mock.group <- X

  res <- run_edgeR_Zinger(SCEdata = scData, condition.name = "mock.group", covariates = NULL,
                    expression.unit = "count", prop.zero = 0.99)
  res$df
})
saveRDS(mock.res.edgeR_Zinger.list, ".../mock_comparison/data10x/mock.res.edgeR_Zinger.list.rds")

mock.res.DESeq2_Zinger.list <- lapply(1:itr, function(i){
  print(i)
  scData        <- scNGP.10x_vehicle[, mock.X[, i]]
  scData$Group <- X
  

  res <- run_DESeq2_Zinger(SCEdata = scData, condition.name = "Group", covariates = NULL,
                          expression.unit = "count", prop.zero = 0.99)
  res$df
})
saveRDS(mock.res.DESeq2_Zinger.list, ".../mock_comparison/data10x/mock.res.DESeq2_Zinger.list.rds")


 
boxplot(list(
PIM = rowMeans(sapply(mock.res.pim.list, function(x){
  p <- x$p.value ; names(p) <- x$ID
  as.numeric(p[rownames(scNGP.10x_vehicle)]<0.05)
}), na.rm=TRUE),
PIM_CoxPH = rowMeans(sapply(mock.res.pim.list2, function(x){
  p <- x$p.value ; names(p) <- x$ID
  as.numeric(p[rownames(scNGP.10x_vehicle)]<0.05)
}), na.rm=TRUE),
MAST = rowMeans(sapply(mock.res.MAST.list, function(x){
  p <- x$pval ; names(p) <- rownames(x)
  as.numeric(p[rownames(scNGP.10x_vehicle)]<0.05)
}), na.rm=TRUE),
DESeq2_Zinger=rowMeans(sapply(mock.res.DESeq2_Zinger.list, function(x){
  p <- x$pvalue ; names(p) <- rownames(x)
  as.numeric(p[rownames(scNGP.10x_vehicle)]<0.05)
}), na.rm=TRUE),
edgeR_Zinger=rowMeans(sapply(mock.res.edgeR_Zinger.list, function(x){
  p <- x$pval; names(p) <- rownames(x)
  as.numeric(p[rownames(scNGP.10x_vehicle)]<0.05)
}), na.rm=TRUE)),
        ylim=c(0, 1))
abline(h=0.05, col=2, lty=2)


hist(mock.res.pim.list[[2]]$p.value)
hist(mock.res.pim.list2[[2]]$p.value)
hist(mock.res.MAST.list[[2]]$pval)
hist(mock.res.edgeR_Zinger.list[[2]]$pval)
hist(mock.res.DESeq2_Zinger.list[[2]]$pvalue)
 
