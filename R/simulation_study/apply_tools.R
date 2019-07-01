source("tools_wrap_functions.R")

# ------------ SPsimSeq simulated data analysis A -------------------------------
sim.data.sp <- readRDS("Simulation_study/SPsimA/sp.sim.dataA.rds")
itr <- length(sim.data.sp)

#PIMseq
resPIMSeq <- lapply(1:itr, function(k){
  sim.data <- sim.data.sp[[k]]
  # keep  <- rowSums(counts(sim.data)>1)>5
  # sim.data <- sim.data[which(keep), ] 
  sim.data$logLS <- log(colData(sim.data)$sim.Lib.Size)
  res <- run_PIMseq(SCEdata = sim.data, condition.name = "Group",
                    covariates = c("logLS"), expression.unit = "counts",  prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.data), by.x="ID", by.y="row.names") 
  #res$all.res[, -1] <- sapply(res$all.res[,-1], function(x) as.numeric(as.character(x)))
  #df <- merge(df, res$all.res, by="ID")
  
  df2 <- res$res$augmented.MP
  df2 <- merge(df2, rowData(sim.data), by.x="ID", by.y="row.names") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(resPIMSeq, "Simulation_study/SPsimA/resPIMSeq.rds")

#MAST
resMAST <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                  expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})
saveRDS(resMAST, "Simulation_study/SPsimA/resMAST.rds")


#SAMSeq
resSAMseq <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.sp[[k]], condition.name = "Group",  
                    covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})
saveRDS(resSAMseq, "Simulation_study/SPsimA/resSAMseq.rds")


#edgeR+Zinger
res.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group", 
                          covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})
saveRDS(res.edgeR_Zinger, "Simulation_study/SPsimA/resZinger_EdgeR.rds")

#DESeq2+Zinger
res.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                           covariates = NULL, expression.unit = "count", 
                           prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})

saveRDS(res.DESeq2_Zinger, "Simulation_study/SPsimA/resZinger_DESeq2.rds")





# ------------ Splat simulated data analysis  A -------------------------------
sim.data.nb <- readRDS("Simulation_study/NBsimA/nb.sim.dataA.rds")
itr <- length(sim.data.nb)
 
#PIMseq
resPIMSeq <- lapply(1:itr, function(k){
  sim.data <- sim.data.nb[[k]]
  sim.data$logLS <- log(colData(sim.data)$ExpLibSize)
  res <- run_PIMseq(SCEdata = sim.data, condition.name = "Group",
                    covariates = c("logLS"), expression.unit = "counts",  prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.data), by.x="ID", by.y="Gene") 
  
  
  df2 <- res$res$augmented.MP
  df2 <- merge(df2, rowData(sim.data), by.x="ID", by.y="row.names") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(resPIMSeq, "Simulation_study/NBsimA/resPIMSeq.rds")


#SAMSeq
resSAMseq <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.nb[[k]], condition.name = "Group",  
                    covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(resSAMseq, "Simulation_study/NBsimA/resSAMseq.rds")

#MAST
resMAST <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                  expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(resMAST, "Simulation_study/NBsimA/resMAST.rds")


#edgeR+Zinger
res.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group", 
                          covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(res.edgeR_Zinger, "Simulation_study/NBsimA/resZinger_EdgeR.rds")

#DESeq2+Zinger
res.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                           expression.unit = "count", prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})

saveRDS(res.DESeq2_Zinger, "Simulation_study/NBsimA/resZinger_DESeq2.rds")




# ------------ SPsimSeq simulated data analysis B -------------------------------
sim.data.sp <- readRDS("Simulation_study/SPsimB/sp.sim.dataB.rds")
itr <- length(sim.data.sp)
 
#PIMseq
resPIMSeq <- lapply(1:itr, function(k){
  sim.data <- sim.data.sp[[k]]
  # keep  <- rowSums(counts(sim.data)>1)>5
  # sim.data <- sim.data[which(keep), ] 
  sim.data$logLS <- log(colData(sim.data)$sim.Lib.Size)
  res <- run_PIMseq(SCEdata = sim.data, condition.name = "Group",
                    covariates = c("logLS"), expression.unit = "counts",  prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.data), by.x="ID", by.y="row.names") 
  #res$all.res[, -1] <- sapply(res$all.res[,-1], function(x) as.numeric(as.character(x)))
  #df <- merge(df, res$all.res, by="ID")
  
  df2 <- res$res$augmented.MP
  df2 <- merge(df2, rowData(sim.data), by.x="ID", by.y="row.names") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(resPIMSeq, "Simulation_study/SPsimB/resPIMSeq.rds")

#MAST
resMAST <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                  expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})
saveRDS(resMAST, "Simulation_study/SPsimB/resMAST.rds")



#edgeR+Zinger
res.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group", 
                          covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})
saveRDS(res.edgeR_Zinger, "Simulation_study/SPsimB/resZinger_EdgeR.rds")

#DESeq2+Zinger
res.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.sp[[k]], condition.name = "Group",
                           covariates = NULL, expression.unit = "count", 
                           prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})

saveRDS(res.DESeq2_Zinger, "Simulation_study/SPsimB/resZinger_DESeq2.rds")

#SAMSeq
resSAMseq <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.sp[[k]], condition.name = "Group",  
                    covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.sp[[k]]), by.x="ID", by.y="row.names") 
  
  list(df=df)
})
saveRDS(resSAMseq, "Simulation_study/SPsimB/resSAMseq.rds")





# ------------ Splat simulated data analysis  B -------------------------------
sim.data.nb <- readRDS("Simulation_study/NBsimB/nb.sim.dataB.rds")
itr <- length(sim.data.nb)

#PIMseq
resPIMSeq <- lapply(1:itr, function(k){
  sim.data <- sim.data.nb[[k]]
  sim.data$logLS <- log(colData(sim.data)$ExpLibSize)
  res <- run_PIMseq(SCEdata = sim.data, condition.name = "Group",
                    covariates = c("logLS"), expression.unit = "counts",  prop.zero = 0.99)  
  df <- merge(res$df, rowData(sim.data), by.x="ID", by.y="Gene") 
  
  
  df2 <- res$res$augmented.MP
  df2 <- merge(df2, rowData(sim.data), by.x="ID", by.y="row.names") 
  
  list(df=df, timing =res$timing, res.augmented= df2)
})
saveRDS(resPIMSeq, "Simulation_study/NBsimB/resPIMSeq.rds")


#SAMSeq
resSAMseq <- lapply(1:itr, function(k){
  res <- run_SAMseq(SCEdata = sim.data.nb[[k]], condition.name = "Group",  
                    covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(resSAMseq, "Simulation_study/NBsimB/resSAMseq.rds")

#MAST
resMAST <- lapply(1:itr, function(k){
  res <- run_MAST(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                  expression.unit = "CPM", covariates = NULL, prop.zero = 0.99)  
  
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(resMAST, "Simulation_study/NBsimB/resMAST.rds")


#edgeR+Zinger
res.edgeR_Zinger <- lapply(1:itr, function(k){
  res <- run_edgeR_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group", 
                          covariates = NULL, expression.unit = "count", prop.zero = 0.99)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})
saveRDS(res.edgeR_Zinger, "Simulation_study/NBsimB/resZinger_EdgeR.rds")

#DESeq2+Zinger
res.DESeq2_Zinger <- lapply(1:itr, function(k){
  res <- run_DESeq2_Zinger(SCEdata = sim.data.nb[[k]], condition.name = "Group",
                           expression.unit = "count", prop.zero = 0.99)
  res$df <- as.data.frame(res$df)
  res$df$ID <- rownames(res$df) 
  df <- merge(res$df, rowData(sim.data.nb[[k]]), by.x="ID", by.y="Gene") 
  
  list(df=df)
})

saveRDS(res.DESeq2_Zinger, "Simulation_study/NBsimB/resZinger_DESeq2.rds")


 
   
