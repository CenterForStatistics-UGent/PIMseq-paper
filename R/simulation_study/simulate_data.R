# --------------  Simulation 1 starting from data A -------------------------
# source data (SMARTer/C1)
library(SingleCellExperiment)
scNGP.polyA <- readRDS("Data/scNGP_polyA.rds") # Data A
dim(scNGP.polyA)
X0 <- colData(scNGP.polyA)$characteristics..treatment
table(X0)

keep <- apply(counts(scNGP.polyA), 1, function(x){
  all(tapply(x==0, X0, sum)>=5)
})
table(keep)

scNGP.polyA2 <- scNGP.polyA[keep,] 
dim(scNGP.polyA2) 
s.data <- scNGP.polyA2
saveRDS(s.data, "Simulation_study/SPsimA/s.dataA.rds")

X <- colData(s.data)$characteristics..treatment
library(SPsimSeq)
itr <- 30

# simulate data
X2 <- ifelse(X=="nutlin",2,1)
sim.data.sp <-  SPsimSeq(n.sim = itr, s.data = s.data, batch = NULL,
                         group = X2, n.genes = 2500, batch.config = 1,
                         group.config = c(0.5, 0.5), tot.samples = 100, pDE = 0.1,
                         model.zero.prob = TRUE, result.format = "SCE", 
                         lfc.thrld = 1,  llStat.thrld = 5, t.thrld = 2.5,  
                         max.frac.zeror.diff = 0.25,  const = 1, seed = 6713)

sim.data.sp[[1]]
saveRDS(sim.data.sp, "Simulation_study/SPsimA/sp.sim.dataA.rds")

#NB simulation from Splatter package
library(splatter)
s.data <- readRDS("Simulation_study/SPsimA/s.dataA.rds")

estimated.params <- splatEstimate(s.data)
set.seed(3125)
seed.set <- sample(3000, itr)
sim.data.nb      <- lapply(1:itr, function(k){
  print(k)
  sim.dat <- splatSimulate(estimated.params, nGenes = 2500, batchCells = 100,
                           group.prob = c(0.5, 0.5), method = "groups", 
                           dropout.type="experiment", de.prob = 0.05,
                           de.facLoc = 2.5, de.facScale = 0.2, batch.facLoc = 0,
                           batch.facScale = 0, seed=seed.set[k])
  rowData(sim.dat)$DE.ind = ifelse(rowData(sim.dat)$DEFacGroup1==1 & rowData(sim.dat)$DEFacGroup2==1, 0, 1)
  sim.dat
})  
saveRDS(estimated.params, "Simulation_study/NBsimA/nb.sim.est.paramsA.rds")
saveRDS(sim.data.nb, "Simulation_study/NBsimA/nb.sim.dataA.rds")


# --------------  Simulation 2 starting from data B -------------------------
# source data (Chromium)
scNGP.10x <- readRDS("Data/NGP_10x_data.rds") # choose data B 
dim(scNGP.10x)
table(colData(scNGP.10x)$treatment)
table(isSpike(scNGP.10x))

X <- colData(scNGP.10x)$treatment
set.seed(6107)
sel.cells <- c(sample(colnames(scNGP.10x[X=="control"]), 500),
               sample(colnames(scNGP.10x[X=="treated"]), 500))
scNGP.10x2 <-scNGP.10x[, sel.cells]
dim(scNGP.10x2) 

keep <- apply(counts(scNGP.10x2), 1, function(x){
  all(tapply(x==0, scNGP.10x2$treatment, sum)>=20)
})
table(keep)
scNGP.10x2 <- scNGP.10x2[keep,]
dim(scNGP.10x2) 

s.data = scNGP.10x2
saveRDS(s.data, "Simulation_study/SPsimB/s.dataB.rds")

X <- colData(s.data)$treatment
X2 <- ifelse(X=="treated",2, 1)
library(SPsimSeq)
itr <- 30
# simulate data
sim.data.sp <-  SPsimSeq(n.sim = itr, s.data = s.data, batch = NULL,
                         group = X2, n.genes = 2500, batch.config = 1,
                         group.config = c(0.5, 0.5), tot.samples = 200, pDE = 0.1,
                         model.zero.prob = TRUE, result.format = "SCE", 
                         lfc.thrld = 0.5,  llStat.thrld = 100, t.thrld = 3.0,  
                         max.frac.zeror.diff = 0.3,  const = 1, seed=5474)

sim.data.sp[[1]]
saveRDS(sim.data.sp, "Simulation_study/SPsimB/sp.sim.dataB.rds")

#NB simulation from Splatter package
library(splatter) 
s.data <- readRDS("Simulation_study/SPsimB/s.dataB.rds")

estimated.params <- splatEstimate(s.data)
set.seed(3125)
seed.set <- sample(3000, itr)
sim.data.nb      <- lapply(1:itr, function(k){
  print(k)
  sim.dat <- splatSimulate(estimated.params, nGenes = 2500, batchCells = 200,
                           group.prob = c(0.5, 0.5), method = "groups", 
                           dropout.type="experiment", de.prob = 0.05,
                           de.facLoc = 1.5, de.facScale = 0.2, batch.facLoc = 0,
                           batch.facScale = 0, seed=seed.set[k])
  rowData(sim.dat)$DE.ind = ifelse(rowData(sim.dat)$DEFacGroup1==1 & rowData(sim.dat)$DEFacGroup2==1, 0, 1)
  sim.dat
})  
saveRDS(estimated.params, "Simulation_study/NBsimB/nb.sim.est.paramsB.rds")
saveRDS(sim.data.nb, "Simulation_study/NBsimB/nb.sim.dataB.rds")


