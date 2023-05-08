load("docs/01-beginner-tmb/exercises/stownet.RData")

Dat <- list(
  Codend = stownet$Codend,
  Cover = stownet$Cover,
  Length = stownet$Lenclass
)

Par <- list(
  beta0 = 0,
  beta1 = 0
)
library(TMB)
compile("docs/01-beginner-tmb/exercises/codend.cpp")
dyn.load(dynlib("docs/01-beginner-tmb/exercises/codend"))

obj <- MakeADFun(Dat, Par, DLL = "codend")
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
opt$objective

sdr <- sdreport(obj)
sdr
summary(sdr, "report")

