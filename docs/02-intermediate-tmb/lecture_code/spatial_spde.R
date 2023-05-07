library(INLA)
library(TMB)
setwd("C:/UAUCK/Conferences/ISEC_2018/TMBWorkshop/SpatialExample")

## -- Simulate data -- ##

#Simulate locations
set.seed(888)
Loc <- matrix(runif(200,0,100),ncol=2)
plot(Loc)

#Generate INLA mesh
mesh <- inla.mesh.create(Loc)
mesh <- inla.mesh.2d(loc=Loc, max.edge = 20, offset = 20)
plot(mesh,asp=1);points(Loc, pch=16)

#Create SPDE object
spde <- inla.spde2.matern(mesh)

# Compile
compile("spatial_spde.cpp" )
dyn.load( dynlib("spatial_spde") )

#Set true parameters
Beta0 <- 1
Alpha <- 5
Kappa <- sqrt(8)/30 #based on range = 30
sig2.omega <- 1
Tau <- sqrt(1/(4*pi*Kappa^2*sig2.omega))

#Set up TMB model
Data <- list(Y_i = rep(0,nrow(Loc)), v_i = mesh$idx$loc-1)
Data$spde <- spde$param.inla[c("M0","M1","M2")]

Params <- list(beta0 = Beta0, ln_alpha = log(Alpha),
              ln_kappa = log(Kappa), ln_tau = log(Tau),
              omega_v = rep(0, mesh$n))

obj <- MakeADFun(Data, Params, DLL = "spatial_spde")
#Simulate data
simdata <- replicate(100, sim <- obj$simulate(obj$par))

set.seed(123)
simdata <- obj$simulate(obj$par)
hist(simdata$Y_i, breaks = 30)

## -- Fit TMB model to estimate parameters -- ##
#Set up TMB model
Data <- list(Y_i = simdata$Y_i, v_i = mesh$idx$loc-1)
Data$spde <- spde$param.inla[c("M0","M1","M2")]
Params <- list(beta0 = 0, ln_alpha = 0,
               ln_kappa = 0, ln_tau = 0,
               omega_v = rep(0, mesh$n))
obj <- MakeADFun(Data, Params, random = "omega_v", DLL = "spatial_spde")
opt <- nlminb(obj$par, obj$fn, obj$gr)

#Report back results
sdr <- sdreport(obj)
c(Beta0, log(Alpha), log(Kappa), log(Tau))
sdr
report <- obj$report()
plot(simdata$Y_i, exp(report$ln_mean))
dyn.unload(dynlib("spatial_spde"))
