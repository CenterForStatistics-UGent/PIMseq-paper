#Single Cell RNA-seq data analysis piplines 

#Preparing input
#source("C:/Users/Alemu/Dropbox/scRNA-seq study/Main scRNA-seq study data analysis/scRNA-research/PIMSeq/main_codes/wrapFunction.R")

prepareInput <- function(SCEdata, condition.name, expression.unit, covariates, ncores, p=0.99)
{
  message("Preparing inputs ...")
  
  #suppressPackageStartupMessages(library(monocle))
  
  if(class(SCEdata) != "SingleCellExperiment"){stop("Input is not a 'SingleCellExperiment' class")}
  
  
  if(!(condition.name %in% colnames(colData(SCEdata)))){
    stop(paste(condition.name, "is not found in colData(SCEdata)"))}
  if(!is.null(covariates) & !all(covariates %in% colnames(colData(SCEdata)))){
    stop(paste(covariates[which(!(covariates %in% colnames(colData(SCEdata))))], 
               "is/are not found in colData(SCEdata)"))
  }
  
  
  #Filtration of low expressed genes
  keep.tag <- apply(counts(SCEdata), 1, function(y){
    mean(y==0) <= p
  })
  
  SCEdata <- SCEdata[rownames(SCEdata) %in% names(keep.tag[keep.tag==TRUE]), ]
  
  
  L <- list()
  
  L[["conditions"]] <- condition.name
  L[["covariates"]] <- covariates
  
  col.info <- data.frame(colData(SCEdata)[, condition.name],  
                         colData(SCEdata)[, covariates], 
                         row.names = rownames(colData(SCEdata)))
  colnames(col.info) <- c(condition.name, covariates)
  L[["col.data"]] <- col.info
  
  L[["n.cores"]] <- ncores
  #L[["total.count"]] <- colSums(counts(SCEdata))
  
  
  if(expression.unit=="count"){
    L[["expr.data"]]       <- counts(SCEdata)
    L[["expression.unit"]] <- expression.unit
    return(L)
  }
  else if(expression.unit=="TPM") {
    if(is.null(rowData(SCEdata)$bp_length)) {stop("Gene length is not provided.")}
    L[["expr.data"]] <- calculateTPM(SCEdata, effective_length = rowData(SCEdata)$bp_length, 
                                     calc_from = "counts")
    L[["expression.unit"]] = expression.unit
    return(L)
  }
  else if(expression.unit=="CPM"){
    L[["expr.data"]] <- edgeR::cpm(counts(SCEdata), normalized.lib.sizes = TRUE)
    L[["expression.unit"]] = expression.unit
    return(L)
  }
  else if(expression.unit=="Census"){
    message("TPM is calculated from raw counts and effective gene length")
    
    if(is.null(rowData(SCEdata)$bp_length)){
      stop("Calculation of Census is failed because gene length is not provided.")}
    
    L[["tpm"]] <- calculateTPM(SCEdata, effective_length = rowData(SCEdata)$bp_length, 
                               calc_from = "counts")
    cds <- newCellDataSet(L$tpm, 
                          phenoData = new("AnnotatedDataFrame", 
                                          data = data.frame(condition = L$col.data[, L$conditions], 
                                                            row.names = colnames(L$tpm))))
    L[["expr.data"]] <- monocle::relative2abs(cds)
    L[["expression.unit"]] = expression.unit
    return(L)
  }
  else{
    stop("Expression unit is not valid. Please choose among 'count', 'TPM', 'CPM', or 'Census'")
  }
}

#PIMSeq
library(PIMseq)
run_PIMseq <- function(SCEdata, condition.name=NULL, covariates=NULL, 
                     expression.unit="counts", prop.zero=0.99, coxph.aprox=FALSE,
                     ncores=1) {
  message("PIMSeq")
  SCEdata <- SCEdata[rowMeans(counts(SCEdata)==0)<prop.zero, ]
  #L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name, p = prop.zero,
  #                  expression.unit = expression.unit, covariates = covariates, ncores = ncores)
  
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      pim.res <- PIMSeq(SCExp = SCEdata, condition = condition.name, nuisance.vars = covariates,
                        assay.name = expression.unit, coxph.aprox=coxph.aprox)
    })
    
    #hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
    
    list(session_info = session_info,
         timing = timing, 
         res = pim.res,
         df = pim.res$test.contrasts)
  }, error = function(e) {
    "PIMseq results could not be calculated"
    list(session_info = session_info)
  })
}


# MAST
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(edgeR))
run_MAST <- function(SCEdata, condition.name=NULL, covariates=NULL,
                     expression.unit, prop.zero=0.99,
                     ncores=NULL) {
  message("MAST, CPM (including detection rate)")

  L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name, p = prop.zero,
                     expression.unit = expression.unit, covariates = covariates, ncores = ncores)

  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      #stopifnot(all(names(L$col.data[,L$conditions]) == colnames(L$expr.data)))
      grp <- L$col.data[,L$conditions]
      names(grp) <- colnames(L$expr.data)
      cdr <- scale(colMeans(L$expr.data > 0))
      batch <- L$col.data[, L$covariates]
      # dge <- DGEList(counts = L$expr.data)
      # dge <- edgeR::calcNormFactors(dge)
      # cpms <- cpm(dge)
      sca <- FromMatrix(exprsArray = log2(L$expr.data + 1),
                        cData = data.frame(wellKey = names(grp),
                                           grp = grp, cdr = cdr, batch=batch))
      if(!is.null(covariates)){
        zlmdata <- zlm(~cdr + batch + grp, sca, silent = TRUE)
      }
      else{
        zlmdata <- zlm(~cdr + grp, sca)
      }
      mast <- lrTest(zlmdata, "grp")
    })

    #hist(mast[, "hurdle", "Pr(>Chisq)"], 50)

    list(session_info = session_info,
         timing = timing,
         res = list(mast=mast, zlmdata=zlmdata),
         df = data.frame(test.stat = mast[, "hurdle", "lambda"],
                         pval = mast[, "hurdle", "Pr(>Chisq)"],
                         qval = p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"])),
         df2 = data.frame(test.stat = mast[, "disc", "lambda"],
                         pval = mast[, "disc", "Pr(>Chisq)"],
                         qval = p.adjust(mast[, "disc", "Pr(>Chisq)"], method = "BH"),
                         row.names = names(mast[, "disc", "Pr(>Chisq)"])))
  }, error = function(e) {
    "MASTcpmDetRate results could not be calculated"
    list(session_info = session_info)
  })
}

# ## Classical wilcoxon test
# suppressPackageStartupMessages(library(edgeR))
# run_Wilcoxon <- function(SCEdata, condition.name=NULL, covariates=NULL,
#                          expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL) {
#   message("Wilcoxon ")
#   
#   if(!is.null(covariates)){message("Covariates are ignored.")}
#   L <- prepare_input(SCEdata = SCEdata, condition.name = condition.name, 
#                      expression.unit = expression.unit, covariates = NULL, ncores = ncores)
#   
#   session_info <- sessionInfo()
#   timing <- system.time({
#     tmm <- edgeR::calcNormFactors(L$expr.data)
#     tpmtmm <- edgeR::cpm(L$expr.data, lib.size = tmm * colSums(L$expr.data))
#     idx <- 1:nrow(tpmtmm)
#     names(idx) <- rownames(tpmtmm)
#     wilcox_p <- sapply(idx, function(i) {
#       wilcox.test(tpmtmm[i, ] ~ L$col.data[, L$conditions])$p.value
#     })
#   })
#   
#   #hist(wilcox_p, 50)
#   
#   list(session_info = session_info,
#        timing = timing,
#        df = data.frame(pval = wilcox_p,
#                        qval = p.adjust(wilcox_p, method = "BH"),
#                        row.names = names(wilcox_p)))
# }

# ## SCDE
# suppressPackageStartupMessages(library(edgeR))
# 
# suppressPackageStartupMessages(library(scde))
# run_SCDE <- function(SCEdata, condition.name=NULL, covariates=NULL,
#                      expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL) {
#   message("scde")
#   L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name,
#                      expression.unit = expression.unit, covariates = NULL, ncores = ncores)
# 
#   session_info <- sessionInfo()
# 
#   timing <- system.time({
#     intcount <- apply(L$expr.data, 2, function(x) {storage.mode(x) <- 'integer'; x})
#     o.ifm <- scde.error.models(counts = intcount, groups = as.factor(unname(L$col.data[, L$conditions])), n.cores = L$n.cores,
#                                save.crossfit.plots = FALSE, save.model.plots = FALSE,
#                                verbose = 0, min.size.entries = min(2000, nrow(L$count) - 1))
#     valid.cells <- o.ifm$corr.a > 0
#     #table(valid.cells)
#     o.ifm <- o.ifm[valid.cells, ]
#     o.prior <- scde.expression.prior(models = o.ifm, counts = intcount[, valid.cells],
#                                      length.out = 400, show.plot = FALSE)
#     grp <- factor(L$col.data[which(valid.cells), L$conditions])
#     names(grp) <- rownames(o.ifm)
#     ediff <- scde.expression.difference(o.ifm, intcount[, valid.cells], o.prior,
#                                         groups = grp, n.randomizations = 100,
#                                         n.cores = L$n.cores, verbose = 0)
#     p.values <- 2*pnorm(abs(ediff$Z), lower.tail = FALSE)
#     p.values.adj <- 2*pnorm(abs(ediff$cZ), lower.tail = FALSE)
#   })
# 
#   #hist(p.values, 50)
#   #hist(p.values.adj, 50)
# 
#   list(session_info = session_info,
#        timing = timing,
#        res = ediff,
#        df = data.frame(pval = p.values,
#                        qval = p.values.adj,
#                        score = abs(ediff$Z),
#                        row.names = rownames(ediff)))
# }


# #scDD
# suppressPackageStartupMessages(library(SummarizedExperiment))
# suppressPackageStartupMessages(library(scran))
# suppressPackageStartupMessages(library(scDD))
# run_scDD <- function(SCEdata, condition.name=NULL, covariates=NULL,
#                       expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL) {
#   message("scDD")
# 
#   L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name,
#                      expression.unit = expression.unit, covariates = covariates, ncores = ncores)
#   if("total_features" %in% colnames(L$col.data)){
#     L$col.data[, "total_features"] <- scale(L$col.data[, "total_features"])
#   }
#   session_info <- sessionInfo()
#   tryCatch({
#     timing <- system.time({
#       # scDatList <- list()
#       # (groups <- unique(L$condt))
#       # for (i in 1:length(groups)) {
#       #   scDatList[[paste0("G", i)]] <- as.matrix(L$count[, which(L$condt == groups[i])])
#       # }
#       scDatList <- SingleCellExperiment(assays=list(counts=L$expr.data),
#                                         colData=L$col.data)
#       datNorm.scran <- scDD::preprocess(scDatList,
#                                         condition = L$conditions,
#                                         zero.thresh = 1, scran_norm = TRUE)
#       condition <- L$col.data[colnames(datNorm.scran), L$conditions]
#       condition <- as.numeric(as.factor(L$col.data[, L$conditions]))
#       names(condition) <- colnames(datNorm.scran)
# 
#       SDSumExp <-SingleCellExperiment(assays = list("normcounts" = normcounts(datNorm.scran)),
#                                        colData = data.frame(condition))
#       prior_param <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
#       scd <- scDD(SDSumExp, prior_param = prior_param, testZeroes = FALSE,
#                   condition = "condition", min.size = 3, min.nonzero = NULL)
#       res <- scDD::results(scd)
#     })
# 
#     #hist(res$nonzero.pvalue, 50)
#     #hist(res$nonzero.pvalue.adj, 50)
# 
#     list(session_info = session_info,
#          timing = timing,
#          res = res,
#          df = data.frame(pval = res$nonzero.pvalue,
#                          qval = res$nonzero.pvalue.adj,
#                          row.names = rownames(res)))
#   }, error = function(e) {
#     "scDD results could not be calculated"
#     list(session_info = session_info)
#   })
# }

#SAMSeq
suppressPackageStartupMessages(library(samr))
run_SAMseq <- function(SCEdata, condition.name=NULL, covariates=NULL, prop.zero=0.95,
                       expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL,
                       resp.type = "Two class unpaired") {
  message("SAMseq")
  message("Any covariate is ignored.")
  L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name, p = prop.zero,
                     expression.unit = expression.unit, covariates = covariates, ncores = ncores)

  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      x=as.numeric(as.factor(L$col.data[, L$conditions]))
      SAMseq.test <- SAMseq(round(L$expr.data), x,
                            resp.type = resp.type,
                            geneid = rownames(L$expr.data), genenames = rownames(L$expr.data),
                            nperms = 100, nresamp = 20, fdr.output = 1)
      SAMseq.result.table <- rbind(SAMseq.test$siggenes.table$genes.up,
                                   SAMseq.test$siggenes.table$genes.lo)
      SAMseq.FDR <- matrix(NA,  nrow=nrow(L$expr.data), ncol=3)#rep(NA, nrow(L$expr.data))
      SAMseq.FDR[match(SAMseq.result.table[, 1], rownames(L$expr.data)), 1] <-
        as.numeric(SAMseq.result.table[, "Score(d)"])
      SAMseq.FDR[match(SAMseq.result.table[, 1], rownames(L$expr.data)), 2] <-
        as.numeric(SAMseq.result.table[, "Fold Change"])
      SAMseq.FDR[match(SAMseq.result.table[, 1], rownames(L$expr.data)), 3] <-
        as.numeric(SAMseq.result.table[, "q-value(%)"])/100
      colnames(SAMseq.FDR) <- c("Score(d)", "FoldChange", "qval")
      rownames(SAMseq.FDR) <- rownames(L$expr.data)
      SAMseq.FDR[is.na(SAMseq.FDR[, 3]), 3] <- 1
      #SAMseq.FDR[SAMseq.FDR[, 2]==0, 3] <- 1
    })

    #hist(SAMseq.FDR, 50)

    list(session_info = session_info,
         timing = timing,
         res = SAMseq.result.table,
         df = as.data.frame(SAMseq.FDR))
  }, error = function(e) {
    "SAMseq results could not be calculated"
    list(session_info = session_info)
  })
}

# ## limma Voom
# suppressPackageStartupMessages(library(limma))
# suppressPackageStartupMessages(library(edgeR))
# 
# run_voomlimma <- function(SCEdata, condition.name=NULL, covariates=NULL,
#                           expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL) {
#   message("voomlimma")
# 
#   L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name,
#                      expression.unit = expression.unit, covariates = covariates, ncores = ncores)
#   if("total_features" %in% colnames(L$col.data)){
#     L$col.data[, "total_features"] <- scale(L$col.data[, "total_features"])
#   }
# 
#   session_info <- sessionInfo()
#   timing <- system.time({
#     dge <- DGEList(L$expr.data, group = L$col.data[, L$conditions])
#     dge <- calcNormFactors(dge)
# 
#     ff = colnames(L$col.data)
#     form.ff = sapply(ff, function(f) {
#       paste0("L$col.data$", f)
#     }, USE.NAMES = FALSE)
#     design <- model.matrix(as.formula(paste("~", paste(form.ff, collapse = "+"))))
# 
#     vm <- voom(dge, design = design, plot = FALSE)
#     fit <- lmFit(vm, design = design)
#     fit <- eBayes(fit)
#     tt <- topTable(fit, n = Inf, adjust.method = "BH", sort.by = "none")
#   })
# 
#   #hist(tt$P.Value, 50)
#   #hist(tt$adj.P.Val, 50)
# 
#   #limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
#   #plotMD(fit)
# 
#   list(session_info = session_info,
#        timing = timing,
#        tt = tt,
#        df = data.frame(pval = tt$P.Value,
#                        qval = tt$adj.P.Val,
#                        row.names = rownames(tt)))
# }
# 
# 
# ## Zinger +edgeR GLM
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(zingeR))

run_edgeR_Zinger <- function(SCEdata, condition.name=NULL, covariates=NULL, prop.zero=0.95,
                         expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL, LS=NULL) {
  message("edgeRLRT")

  L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name, p=prop.zero,
                     expression.unit = expression.unit, covariates = covariates, ncores = ncores)
  if("total_features" %in% colnames(L$col.data)){
    L$col.data[, "total_features"] <- scale(L$col.data[, "total_features"])
  }
  session_info <- sessionInfo()
  timing <- system.time({
    dge <- DGEList(L$expr.data, group = L$col.data[, L$conditions])
    if(!is.null(LS)){
      dge$samples$lib.size <- LS
    }
    dge <- suppressWarnings(calcNormFactors(dge))

    ff = 1:length(colnames(L$col.data))
    form.ff = sapply(ff, function(f) {
      paste("L$col.data[,", f,"]")
    }, USE.NAMES = FALSE)
    design <- model.matrix(as.formula(paste("~", paste(form.ff, collapse = "+"))))

    weights <- zeroWeightsLS(counts=dge$counts, design=design, maxit=300, normalization="TMM") #Zinger
    dge$weights <- weights
    dge=estimateDisp(dge, design)
    fit <- glmFit(dge, design = design)
    lrt <- glmWeightedF(fit, coef=2, independentFiltering = TRUE)
    tt  <- topTags(lrt, n = Inf, sort.by = "none")
  })

  # plotBCV(dge)
  # hist(tt$table$PValue, 50)
  # hist(tt$table$FDR, 50)
  # limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
  # plotSmear(lrt)

  list(session_info = session_info,
       timing = timing,
       tt = tt,
       df = data.frame(pval = tt$table$PValue,
                       qval = tt$table$FDR,
                       row.names = rownames(tt$table)))
}
# # 
# # ## DESeq2 + Zinger
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(zingeR))

run_DESeq2_Zinger <- function(SCEdata, condition.name=NULL, covariates=NULL, prop.zero=0.95,
                       expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL) {
  message("DESeq2")
  L <- prepareInput(SCEdata = SCEdata, condition.name = condition.name, p=prop.zero,
                     expression.unit = expression.unit, covariates = covariates, ncores = ncores)
  if("total_features" %in% colnames(L$col.data)){
    L$col.data[, "total_features"] <- scale(L$col.data[, "total_features"]) 
  }
  #L$col.data[, "Group"] <- as.factor(L$col.data[, "Group"])
  L$col.data[, condition.name] <- as.factor(L$col.data[, condition.name])
  session_info <- sessionInfo()
  timing <- system.time({

      ff = colnames(L$col.data)
      form.ff = sapply(ff, function(f) {
        paste(f)
      }, USE.NAMES = FALSE)

      colData=L$col.data

      design <- as.formula(paste("~", paste(form.ff, collapse = "+")))

      dds <- DESeqDataSetFromMatrix(countData = round(L$expr.data),
                                    colData = L$col.data,
                                    design = design)
      weights <- zeroWeightsLS(counts=L$expr.data, design=model.matrix(design, data=L$col.data), maxit=300,
                               normalization="DESeq2_poscounts", colData=colData,
                               designFormula=design)
      assays(dds)[["weights"]]= weights

      dse = DESeq2::estimateSizeFactors(dds, type="poscounts")
      dse = estimateDispersions(dse)
      dse = nbinomWaldTest(dse, betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-2)
      res = DESeq2::results(dse, contrast = c(L$conditions, levels(factor(L$col.data[, L$conditions]))[1],
                                       levels(factor(L$col.data[, L$conditions]))[2]),
                            independentFiltering=FALSE)
      #res$qval <- p.adjust(res$pvalue, method = "BH")
      res$padj[is.na(res$padj)] <- 1
    })

    #plotDispEsts(dds)
    #plotMA(res)
    #summary(res)

    list(session_info = session_info,
         timing = timing,
         res = res,
         df = res)
}


#ROTSvoom
# suppressPackageStartupMessages(library(ROTS))
# suppressPackageStartupMessages(library(edgeR))
# 
# run_ROTSvoom <- function(SCEdata, condition.name=NULL, covariates=NULL,
#                          expression.unit= c("count", "TPM", "CPM", "Census"), ncores=NULL) {
#   message("ROTS, voom")
#   L <- prepare_input(SCEdata = SCEdata, condition.name = condition.name, 
#                      expression.unit = expression.unit, covariates = covariates, ncores = ncores)
#   
#   session_info <- sessionInfo()
#   timing <- system.time({
#     stopifnot(all(rownames(L$col.data) == colnames(L$expr.data)))
#     grp <- as.factor(L$col.data[, L$conditions])
#     dge <- DGEList(counts = L$expr.data)
#     dge <- edgeR::calcNormFactors(dge)
#     vm <- voom(dge, design = model.matrix(~ grp))
#     rots <- ROTS(data = vm$E, groups = as.numeric(grp), B = 1000, K = 1000, log = TRUE, seed = 123,
#                  progress =TRUE)
#   })
#   
#   # hist(rots$pvalue, 50)
#   # hist(rots$FDR, 50)
#   # hist(rots$logfc, 50)
#   # print(rots$R)
#   # print(rots$Z)
#   # print(rots$k)
#   # print(rots$a1)
#   # print(rots$a2)
#   
#   list(session_info = session_info,
#        timing = timing,
#        res = rots,
#        df = data.frame(pval = rots$pvalue,
#                        qval = p.adjust(rots$pvalue, method = "BH"),
#                        row.names = rownames(rots$data)))
# }
#------------------------------------------------------------------------

























