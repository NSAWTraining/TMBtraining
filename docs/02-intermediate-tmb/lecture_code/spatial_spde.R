#install INLA from website (not on CRAN or github)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages("inlabru")
library(INLA)
library(TMB)
library(inlabru)
library(ggplot2)

## -- Simulate data -- ##

#Simulate locations
set.seed(888)
Loc <- data.frame(x = runif(100,0,100), 
                  y = runif(100,0,100))
p <- ggplot(Loc, aes(x=x, y=y)) + geom_point() +
  theme_classic()
p

#Generate INLA mesh
mesh <- inla.mesh.create(as.matrix(Loc))
mesh <- inla.mesh.2d(loc=as.matrix(Loc), max.edge = 20, offset = 20)
p + gg(mesh)

#Create SPDE object
spde <- inla.spde2.matern(mesh)

# Compile
compile("docs/02-intermediate-tmb/lecture_code/spatial_spde.cpp" )
dyn.load( dynlib("docs/02-intermediate-tmb/lecture_code/spatial_spde") )

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

#Map estimated RE field
omega.est <- obj$env$parList()$omega_v
proj.grid <- expand.grid(x = 0:100, y = 0:100)
P <- inla.mesh.project(mesh, as.matrix(proj.grid))
proj.grid$z <- as.vector(P$A %*% omega.est)
ggplot(proj.grid, aes(x=x, y=y, color = z)) + geom_point()

