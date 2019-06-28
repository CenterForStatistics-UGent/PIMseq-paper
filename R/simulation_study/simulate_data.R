# --------------  Simulation 1 starting from data A -------------------------
# source data (SMARTer/C1)
scNGP.polyA <- readRDS(choose.files()) # Data A
dim(scNGP.polyA)
X0 <- colData(scNGP.polyA)$characteristics..treatment
table(X0)
Ux <- as.character(unique(X0))


scNGP.polyA2 <- scNGP.polyA[rowSums(counts(scNGP.polyA[, X0==Ux[1]])>1)>=5 & rowSums(counts(scNGP.polyA[, X0==Ux[2]])>1)>=5,] 
dim(scNGP.polyA2) 
s.data <- scNGP.polyA2
saveRDS(".../s.dataA.rds")

X <- colData(s.data)$characteristics..treatment
library(SPsimSeq)
itr <- 30
# simulate data
sim.data.sp <-  SPsimSeq(n.sim = itr, s.data = s.data, batch = NULL,
                         group = X, n.genes = 2500, batch.config = 1,
                         group.config = c(0.5, 0.5), tot.samples = 100, pDE = 0.1,
                         model.zero.prob = TRUE, result.format = "SCE", 
                         lfc.thrld = 1,  llStat.thrld = 5, t.thrld = 2.5,  
                         max.frac.zeror.diff = 0.2,  const = 1)

sim.data.sp[[1]]
saveRDS(sim.data.sp, ".../sp.sim.dataA.rds")

#NB simulation from Splatter package
library(splatter)
s.data <- readRDS(".../s.dataA.rds")

estimated.params <- splatEstimate(s.data)
set.seed(3125)
seed.set <- sample(3000, itr)
sim.data.nb      <- lapply(1:itr, function(k){
  print(k)
  sim.dat <- splatSimulate(estimated.params, nGenes = 2500, batchCells = 100,
                group.prob = c(0.5, 0.5), method = "groups", 
                dropout.type="experiment", de.prob = 0.05,
                de.facLoc = 1, de.facScale = 0.2, batch.facLoc = 0,
                batch.facScale = 0, seed=seed.set[k])
  rowData(sim.dat)$DE.ind = ifelse(rowData(sim.dat)$DEFacGroup1==1 & rowData(sim.dat)$DEFacGroup2==1, 0, 1)
  sim.dat
})  
saveRDS(estimated.params, ".../nb.sim.est.paramsA.rds")
saveRDS(sim.data.nb, ".../nb.sim.dataA.rds")


# --------------  Simulation 2 starting from data B -------------------------
# source data (Chromium)
scNGP.10x <- readRDS(choose.files()) # choose data B 
dim(scNGP.10x)
table(colData(scNGP.10x)$treatment)
table(isSpike(scNGP.10x))

X <- colData(scNGP.10x)$treatment
set.seed(6107)
sel.cells <- c(sample(colnames(scNGP.10x[X=="control"]), 500),
               sample(colnames(scNGP.10x[X=="treated"]), 500))
scNGP.10x2 <-scNGP.10x[, sel.cells]
dim(scNGP.10x2) 

scNGP.10x2 <- scNGP.10x2[(rowSums(counts(scNGP.10x2[, colData(scNGP.10x2)$treatment=="control"])>0)>=20 &
                         rowSums(counts(scNGP.10x2[, colData(scNGP.10x2)$treatment=="treated"])>0)>=20),]
dim(scNGP.10x2) 

s.data = scNGP.10x2
saveRDS(s.data, ".../s.dataB.rds")

X <- colData(s.data)$treatment
library(SPsimSeq)
itr <- 30
# simulate data
sim.data.sp <-  SPsimSeq(n.sim = itr, s.data = s.data, batch = NULL,
                         group = X, n.genes = 2500, batch.config = 1,
                         group.config = c(0.5, 0.5), tot.samples = 200, pDE = 0.1,
                         model.zero.prob = TRUE, result.format = "SCE", 
                         lfc.thrld = 0.5,  llStat.thrld = 100, t.thrld = 3.0,  
                         max.frac.zeror.diff = 0.3,  const = 1)

sim.data.sp[[1]]
saveRDS(sim.data.sp, ".../sp.sim.dataB.rds")

#NB simulation from Splatter package
library(splatter) 
s.data <- readRDS(".../s.dataB.rds")

estimated.params <- splatEstimate(s.data)
set.seed(3125)
seed.set <- sample(3000, itr)
sim.data.nb      <- lapply(1:itr, function(k){
  print(k)
  sim.dat <- splatSimulate(estimated.params, nGenes = 2500, batchCells = 200,
                group.prob = c(0.5, 0.5), method = "groups", 
                dropout.type="experiment", de.prob = 0.05,
                de.facLoc = 1, de.facScale = 0.2, batch.facLoc = 0,
                batch.facScale = 0, seed=seed.set[k])
  rowData(sim.dat)$DE.ind = ifelse(rowData(sim.dat)$DEFacGroup1==1 & rowData(sim.dat)$DEFacGroup2==1, 0, 1)
  sim.dat
})  
saveRDS(estimated.params, ".../nb.sim.est.paramsB.rds")
saveRDS(sim.data.nb, ".../nb.sim.dataB.rds")
