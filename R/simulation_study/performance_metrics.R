# A function to calculate performance metrices
calcFDP_TPR <- function(res, test="qval", t=0, results_in="df"){
  fdr <- seq(0, 0.5, 0.0025)
  FDP_TPR <- as.data.frame(t(sapply(fdr, function(a){
    t <- apply(sapply(res, function(x){
      df    <- as.data.frame(x[[results_in]])
      df$DE <- ifelse(df$DE.ind==1, 1, 0)
      
      if(sum(df[, test]<a)>0){
        fdp <- sum(df[, test]<a & df$DE==0)/sum(df[, test]<a)
        tpp <- sum(df[, test]<a & df$DE==1)/sum(df$DE==1)
        prc <- sum(df[, test]<a & df$DE==1)/(sum(df[, test]<a & df$DE==1) +
                                               sum(df[, test]<a & df$DE==0))
      }
      else{
        fdp <- 0
        tpp <- 0
        prc <- 0
      }
      c(fdp=fdp, tpp=tpp, prc=prc)
    }), 1, mean, trim=t)
  })))
  FDP_TPR$nom.fdr <- fdr
  FDP_TPR
} 
# --------- SPsimSeq simulation A ------------------------------
runPIMSeq  <- readRDS("Simulation_study/SPsimA/resPIMSeq.rds")
runMAST      <- readRDS("Simulation_study/SPsimA/resMAST.rds")
runSAMseq    <- readRDS("Simulation_study/SPsimA/resSAMSeq.rds")
runZinger_edgeR    <- readRDS("Simulation_study/SPsimA/resZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS("Simulation_study/SPsimA/resZinger_DESeq2.rds")

perf.list <- list(
  PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
  #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
  MAST   = calcFDP_TPR(runMAST),
  SAMSeq = calcFDP_TPR(runSAMseq),
  edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
  DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
) 
saveRDS(perf.list, "Simulation_study/SPsimA/perf.list.rds")

# --------- Splat simulation A ------------------------------
runPIMSeq  <- readRDS("Simulation_study/NBsimA/resPIMSeq.rds")
runMAST      <- readRDS("Simulation_study/NBsimA/resMAST.rds")
runSAMseq    <- readRDS("Simulation_study/NBsimA/resSAMSeq.rds")
runZinger_edgeR    <- readRDS("Simulation_study/NBsimA/resZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS("Simulation_study/NBsimA/resZinger_DESeq2.rds")

perf.list <- list(
  PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
  #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
  MAST   = calcFDP_TPR(runMAST),
  SAMSeq = calcFDP_TPR(runSAMseq),
  edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
  DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
) 
saveRDS(perf.list, "Simulation_study/NBsimA/perf.list.rds")

# ----------  Figure 2
perf.list1 <- readRDS("Simulation_study/SPsimA/perf.list.rds") 
names(perf.list1) <- c("PIM", "MAST", "SAMSeq", "edgeR + Zinger","DESeq2 +Zinger")

par(mfrow=c(1,2))
plot(c(0, 1), c(0, 1), type="n", xlab="FDR", ylab="TPR", las=1, xlim=c(0, 0.4), ylim=c(0, 1),
     main="SPsimSeq simulation")
for(i in 1:length(perf.list1)){
  lines(perf.list1[[i]]$fdp, perf.list1[[i]]$tpp, col=i, lwd=3)
  points(perf.list1[[i]]$fdp[perf.list1[[i]]$nom.fdr==0.05], 
         perf.list1[[i]]$tpp[perf.list1[[i]]$nom.fdr==0.05], col=i, pch=19, cex=2)
}
abline(0, 1, lty=3) ; abline(v=0.05, lty=3)
legend("bottomright", names(perf.list1), col=1:length(perf.list1), lty=1, lwd = 2, pch=19)


perf.list1 <- readRDS("Simulation_study/NBsimA/perf.list.rds")  
names(perf.list1) <- c("PIM", "MAST", "SAMSeq", "edgeR + Zinger","DESeq2 +Zinger")

plot(c(0, 1), c(0, 1), type="n", xlab="FDR", ylab="TPR", las=1, xlim=c(0, 0.4), ylim=c(0, 1),
     main="Splat simulation")
for(i in 1:length(perf.list1)){
  lines(perf.list1[[i]]$fdp, perf.list1[[i]]$tpp, col=i, lwd=3)
  points(perf.list1[[i]]$fdp[perf.list1[[i]]$nom.fdr==0.05], 
         perf.list1[[i]]$tpp[perf.list1[[i]]$nom.fdr==0.05], col=i, pch=19, cex=2)
}
abline(0, 1, lty=3) ; abline(v=0.05, lty=3)
legend("bottomright", names(perf.list1), col=1:length(perf.list1), lty=1, lwd = 2, pch=19)





# --------- SPsimSeq simulation B ------------------------------
runPIMSeq  <- readRDS("Simulation_study/SPsimB/resPIMSeq.rds")
runMAST      <- readRDS("Simulation_study/SPsimB/resMAST.rds")
runSAMseq    <- readRDS("Simulation_study/SPsimB/resSAMSeq.rds")
runZinger_edgeR    <- readRDS("Simulation_study/SPsimB/resZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS("Simulation_study/SPsimB/resZinger_DESeq2.rds")

perf.list <- list(
  PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
  #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
  MAST   = calcFDP_TPR(runMAST),
  SAMSeq = calcFDP_TPR(runSAMseq),
  edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
  DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
) 
saveRDS(perf.list, "Simulation_study/SPsimB/perf.list.rds")



# --------- Splat simulation B ------------------------------
runPIMSeq  <- readRDS("Simulation_study/NBsimB/resPIMSeq.rds")
runMAST      <- readRDS("Simulation_study/NBsimB/resMAST.rds")
runSAMseq    <- readRDS("Simulation_study/NBsimB/resSAMSeq.rds")
runZinger_edgeR    <- readRDS("Simulation_study/NBsimB/resZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS("Simulation_study/NBsimB/resZinger_DESeq2.rds")


perf.list <- list(
  PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
  #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
  MAST   = calcFDP_TPR(runMAST),
  SAMSeq = calcFDP_TPR(runSAMseq),
  edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
  DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
) 
saveRDS(perf.list, "Simulation_study/NBsimB/perf.list.rds")

# ----------  Figure 3
perf.list1 <- readRDS("Simulation_study/SPsimB/perf.list.rds") 
names(perf.list1) <- c("PIM", "MAST", "SAMSeq", "edgeR + Zinger","DESeq2 +Zinger")

par(mfrow=c(1,2))
plot(c(0, 1), c(0, 1), type="n", xlab="FDR", ylab="TPR", las=1, xlim=c(0, 0.4), ylim=c(0, 1),
     main="Semi-parametric simulation")
for(i in 1:length(perf.list1)){
  lines(perf.list1[[i]]$fdp, perf.list1[[i]]$tpp, col=i, lwd=3)
  points(perf.list1[[i]]$fdp[perf.list1[[i]]$nom.fdr==0.05], 
         perf.list1[[i]]$tpp[perf.list1[[i]]$nom.fdr==0.05], col=i, pch=19, cex=2)
}
abline(0, 1, lty=3) ; abline(v=0.05, lty=3)
legend("bottomright", names(perf.list1), col=1:length(perf.list1), lty=1, lwd = 2, pch=19)


perf.list1 <- readRDS("Simulation_study/NBsimB/perf.list.rds")  
names(perf.list1) <- c("PIM", "MAST", "SAMSeq", "edgeR + Zinger","DESeq2 +Zinger")

plot(c(0, 1), c(0, 1), type="n", xlab="FDR", ylab="TPR", las=1, xlim=c(0, 0.4), ylim=c(0, 1),
     main="Negative binomial simulation")
for(i in 1:length(perf.list1)){
  lines(perf.list1[[i]]$fdp, perf.list1[[i]]$tpp, col=i, lwd=3)
  points(perf.list1[[i]]$fdp[perf.list1[[i]]$nom.fdr==0.05], 
         perf.list1[[i]]$tpp[perf.list1[[i]]$nom.fdr==0.05], col=i, pch=19, cex=2)
}
abline(0, 1, lty=3) ; abline(v=0.05, lty=3)
legend("bottomright", names(perf.list1), col=1:length(perf.list1), lty=1, lwd = 2, pch=19)






