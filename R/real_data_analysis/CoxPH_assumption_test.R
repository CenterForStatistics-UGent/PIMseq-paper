
# -------------- C1 data ----------------
scNGP.C1_dat <- readRDS("Data/DataA.rds")
scNGP.C1_filtered <- scNGP.C1_dat$filtered

library(survival)

A <- as.factor(scNGP.C1_filtered$treatment)
x <- log(colSums(counts(scNGP.C1_filtered)))

p.val.PH <- t(sapply(1:nrow(scNGP.C1_filtered), function(i){
  print(i)
  y <- as.numeric(counts(scNGP.C1_filtered)[i,])
  status <- rep(1,length(y))
  res.cox <- coxph(Surv(y, status)~A+x, ties="efron")
  test.ph <- cox.zph(res.cox)
  test.ph$table[,3]
}))
tail(p.val.PH)
par(mfrow=c(1,3))
hist(p.val.PH[,1], main="A")
hist(p.val.PH[,2], main="X")
hist(p.val.PH[,3], main="Global")

library(fdrtool)
frac.H0 <- t(apply(p.val.PH, 2, function(p){
  est <- fdrtool(p, statistic="pvalue", cutoff.method="fndr", plot=FALSE)
  est$param[1, c("eta0", "eta0.SE")]
}))
frac.H0


# ---------------- Chromium data --------------------
library(SingleCellExperiment)
scNGP.10x_dat <- readRDS("Data/DataB.rds")
scNGP.10x_filtered <- scNGP.10x_dat$filtered


library(survival)

A <- as.factor(scNGP.10x_filtered$treatment)
x <- log(colSums(counts(scNGP.10x_filtered)))

p.val.PH <- t(sapply(1:nrow(scNGP.10x_filtered), function(i){
  print(i)
  y <- as.numeric(counts(scNGP.10x_filtered)[i,])
  status <- rep(1,length(y))
  res.cox <- coxph(Surv(y, status)~A+x, ties="efron")
  test.ph <- cox.zph(res.cox)
  test.ph$table[,3]
}))
tail(p.val.PH)
par(mfrow=c(1,3))
hist(p.val.PH[,1], main="A")
hist(p.val.PH[,2], main="X")
hist(p.val.PH[,3], main="Global")

library(fdrtool)
frac.H0 <- t(apply(p.val.PH, 2, function(p){
  est <- fdrtool(p, statistic="pvalue", cutoff.method="fndr", plot=FALSE)
  est$param[1, c("eta0", "eta0.SE")]
}))
frac.H0






