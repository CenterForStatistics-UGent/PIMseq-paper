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
runPIMSeq  <- readRDS(".../SPsimA/runPIMSeq.rds")
runMAST      <- readRDS(".../SPsimA/runMAST.rds")
runSAMseq    <- readRDS(".../SPsimA/runSAMSeq.rds")
runZinger_edgeR    <- readRDS(".../SPsimA/runZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS(".../SPsimA/runZinger_DESeq2.rds")

perf.list <- list(
      PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
      #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
      MAST   = calcFDP_TPR(runMAST),
      SAMSeq = calcFDP_TPR(runSAMseq),
      edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
      DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
   ) 
saveRDS(perf.list, "/SPsimA/perf.list.rds")

# --------- Splat simulation A ------------------------------
runPIMSeq  <- readRDS(".../NBsimA/runPIMSeq.rds")
runMAST      <- readRDS(".../NBsimA/runMAST.rds")
runSAMseq    <- readRDS(".../NBsimA/runSAMSeq.rds")
runZinger_edgeR    <- readRDS(".../NBsimA/runZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS(".../NBsimA/runZinger_DESeq2.rds")

perf.list <- list(
      PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
      #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
      MAST   = calcFDP_TPR(runMAST),
      SAMSeq = calcFDP_TPR(runSAMseq),
      edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
      DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
   ) 
saveRDS(perf.list, "/NBsimA/perf.list.rds")



# --------- SPsimSeq simulation B ------------------------------
runPIMSeq  <- readRDS(".../SPsimB/runPIMSeq.rds")
runMAST      <- readRDS(".../SPsimB/runMAST.rds")
runSAMseq    <- readRDS(".../SPsimB/runSAMSeq.rds")
runZinger_edgeR    <- readRDS(".../SPsimB/runZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS(".../SPsimB/runZinger_DESeq2.rds")

perf.list <- list(
      PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
      #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
      MAST   = calcFDP_TPR(runMAST),
      SAMSeq = calcFDP_TPR(runSAMseq),
      edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
      DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
   ) 
saveRDS(perf.list, "/SPsimB/perf.list.rds")



# --------- Splat simulation B ------------------------------
runPIMSeq  <- readRDS(".../NBsimB/runPIMSeq.rds")
runMAST      <- readRDS(".../NBsimB/runMAST.rds")
runSAMseq    <- readRDS(".../NBsimB/runSAMSeq.rds")
runZinger_edgeR    <- readRDS(".../NBsimB/runZinger_EdgeR.rds")
runZinger_DESeq2   <- readRDS(".../NBsimB/runZinger_DESeq2.rds")

perf.list <- list(
      PIMSeq1A = calcFDP_TPR(runPIMSeq, test = "p.adjusted"),
      #PIMSeq1A_marg    = calcFDP_TPR(runPIMSeq, test = "p.adjusted", results_in="res.augmented"), 
      MAST   = calcFDP_TPR(runMAST),
      SAMSeq = calcFDP_TPR(runSAMseq),
      edgeR_Zinger  = calcFDP_TPR(runZinger_edgeR),
      DESeq2_Zinger = calcFDP_TPR(runZinger_DESeq2, test = "padj")
   ) 
saveRDS(perf.list, "/NBsimB/perf.list.rds")
