
source(".../tools_wrap_functions.R")

# ------------ SPsimSeq simulated data analysis A -------------------------------
sim.data.sp <- readRDS(".../sp.sim.dataA.rds")
itr <- length(sim.data.sp)

#PIMseq
runPIMSeq1_A <- lapply(1:itr, function(k){
  sim.data <- sim.data.sp[[k]]
  # keep  <- rowSums(counts(sim.data)>1)>5
  # sim.data <- sim.data[which(keep), ]
  #plot(colSums(counts(sim.data)), colData(sim.data)$simLibeSize)
  sim.dat$logLS <- log(colData(sim.dat)$ExpLibSize)
  res <- run_PIMSeq(SCEdata = sim.data, condition.name = "Group",
                     covariates = c("logLS"), expression.unit = "count",  prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.data), by.x="ID", by.y="GeneID") 
  res$all.res[, -1] <- sapply(res$all.res[,-1], function(x) as.numeric(as.character(x)))
  df <- merge(df, res$all.res, by="ID")
  
  df2 <- res$res.augmented 
  df2 <- merge(df2, rowData(sim.data), by.x="ID", by.y="GeneID") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(runPIMSeq1_A, ".../SPsimA/runPIMSeq1.rds")

#MAST
runMAST_A <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                     expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 
  
  list(df=df)
})
saveRDS(runMAST_A, ".../SPsimA/runMAST.rds")


#SAMSeq
runSAMseq_A <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.sp[[k]], condition.name = "Group",  
                                  covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 

  list(df=df)
})
saveRDS(runSAMseq_A, ".../SPsimA/runSAMseq.rds")


#edgeR+Zinger
result.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group", 
                          covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 

  list(df=df)
})
saveRDS(result.edgeR_Zinger, ".../SPsimA/runZinger_EdgeR.rds")

#DESeq2+Zinger
result.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                           covariates = NULL, expression.unit = "count", 
                           prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 
  
  list(df=df)
})

saveRDS(result.DESeq2_Zinger, ".../SPsimA/runZinger_DESeq2.rds")




  
# ------------ Splat simulated data analysis  A -------------------------------
sim.data.nb <- readRDS(".../nb.sim.dataA.rds")
itr <- length(sim.data.nb)

#PIMseq
runPIMSeq1_A <- lapply(1:itr, function(k){
  sim.dat <- sim.data.nb[[k]]
  sim.dat$logLS <- log(colData(sim.dat)$ExpLibSize)
  res <- run_PIMseq(SCEdata = sim.dat, condition.name = "Group", 
                    covariates = c("logLS"), prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.dat), by.x="ID", by.y="Gene")  
  df2 <- res$res$augmented.MP
  df2 <- merge(df2, rowData(sim.dat), by.x="ID", by.y="Gene") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(runPIMSeq1_A, ".../NBsimA/runPIMSeq1.rds")

#SAMSeq
runSAMseq_A <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.nb[[k]], condition.name = "Group",  
                                  covariates = NULL, expression.unit = "count", prop.zero = 0.95)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 

  list(df=df)
})
saveRDS(runSAMseq_A, ".../NBsimA/runSAMseq.rds")

#MAST
runMAST_A <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                     expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(runMAST_A, ".../NBsimA/runMAST.rds")


#edgeR+Zinger
result.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group", 
                                        covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 

  list(df=df)
})
saveRDS(result.edgeR_Zinger, ".../NBsimA/runZinger_EdgeR.rds")

#DESeq2+Zinger
result.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                           expression.unit = "count", prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})

saveRDS(result.DESeq2_Zinger, ".../NBsimA/runZinger_DESeq2.rds")

  
  
  
  # ------------ SPsimSeq simulated data analysis A -------------------------------
sim.data.sp <- readRDS(".../sp.sim.dataB.rds")
itr <- length(sim.data.sp)

#PIMseq
runPIMSeq1_A <- lapply(1:itr, function(k){
  sim.data <- sim.data.sp[[k]]
  # keep  <- rowSums(counts(sim.data)>1)>5
  # sim.data <- sim.data[which(keep), ]
  #plot(colSums(counts(sim.data)), colData(sim.data)$simLibeSize)
  sim.dat$logLS <- log(colData(sim.dat)$ExpLibSize)
  res <- run_PIMSeq(SCEdata = sim.data, condition.name = "Group",
                     covariates = c("logLS"), expression.unit = "count",  prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.data), by.x="ID", by.y="GeneID") 
  res$all.res[, -1] <- sapply(res$all.res[,-1], function(x) as.numeric(as.character(x)))
  df <- merge(df, res$all.res, by="ID")
  
  df2 <- res$res.augmented 
  df2 <- merge(df2, rowData(sim.data), by.x="ID", by.y="GeneID") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(runPIMSeq1_A, ".../SPsimB/runPIMSeq.rds")

#MAST
runMAST_A <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                     expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 
  
  list(df=df)
})
saveRDS(runMAST_A, ".../SPsimB/runMAST.rds")


#SAMSeq
runSAMseq_A <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.sp[[k]], condition.name = "Group",  
                                  covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 

  list(df=df)
})
saveRDS(runSAMseq_A, ".../SPsimB/runSAMseq.rds")


#edgeR+Zinger
result.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group", 
                          covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 

  list(df=df)
})
saveRDS(result.edgeR_Zinger, ".../SPsimA/runZinger_EdgeR.rds")

#DESeq2+Zinger
result.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                           covariates = NULL, expression.unit = "count", 
                           prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="GeneID") 
  
  list(df=df)
})

saveRDS(result.DESeq2_Zinger, ".../SPsimB/runZinger_DESeq2.rds")




  
# ------------ Splat simulated data analysis  B -------------------------------
sim.data.nb <- readRDS(".../nb.sim.dataB.rds")
itr <- length(sim.data.nb)

#PIMseq
runPIMSeq1_A <- lapply(1:itr, function(k){
  sim.dat <- sim.data.nb[[k]]
  sim.dat$logLS <- log(colData(sim.dat)$ExpLibSize)
  res <- run_PIMseq(SCEdata = sim.dat, condition.name = "Group", 
                    covariates = c("logLS"), prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.dat), by.x="ID", by.y="Gene")  
  df2 <- res$res$augmented.MP
  df2 <- merge(df2, rowData(sim.dat), by.x="ID", by.y="Gene") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(runPIMSeq1_A, ".../NBsimB/runPIMSeq.rds")

#SAMSeq
runSAMseq_A <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.nb[[k]], condition.name = "Group",  
                                  covariates = NULL, expression.unit = "count", prop.zero = 0.95)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 

  list(df=df)
})
saveRDS(runSAMseq_A, ".../NBsimB/runSAMseq.rds")

#MAST
runMAST_A <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                     expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(runMAST_A, ".../NBsimB/runMAST.rds")


#edgeR+Zinger
result.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group", 
                                        covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 

  list(df=df)
})
saveRDS(result.edgeR_Zinger, ".../NBsimB/runZinger_EdgeR.rds")

#DESeq2+Zinger
result.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                           expression.unit = "count", prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})

saveRDS(result.DESeq2_Zinger, ".../NBsimB/runZinger_DESeq2.rds")


 
   
