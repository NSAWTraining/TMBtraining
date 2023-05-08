#install INLA from website (not on CRAN or github)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#install.packages("inlabru", "mvtnorm", "ggplot2")
library(INLA)
library(TMB)
library(inlabru)
library(ggplot2)
library(mvtnorm)

## -- Simulate data -- ##

## functions for simulating spatial data
cMatern <- function(H, Nu, Kap) {
  ifelse(H > 0, besselK(H*Kap, Nu) * (H*Kap)^Nu, 1) / gamma(Nu) * 2^(1-Nu)
}

# Simulate spatial field
sim_omega <- function(Range, sig2, Dmat){
  Kappa <- sqrt(8)/Range
  n <- dim(Dmat)[1]
  Sigma <- sig2 * cMatern(Dmat,1,Kappa)
  omega <- t(mvtnorm::rmvnorm(1, rep(0,nrow(Sigma)), 
                              sigma = Sigma, method = 'chol'))

  return(omega)
}

#Set up projection grid
proj.grid <- expand.grid(x = 0:100, y = 0:100)


#Set true parameters
Beta0 <- 2
Kappa <- sqrt(8)/50 #based on range = 30
sig2.omega <- 2
Tau <- sqrt(1/(4*pi*Kappa^2*sig2.omega))

#Simulate random spatial field
# set.seed(88)
# omega.sim <- sim_omega(Range = 50, sig2 = 2,
#                        Dmat = as.matrix(dist(as.matrix(proj.grid))))
# spatial.df <- data.frame(x = proj.grid[,1],
#                          y = proj.grid[,2],
#                          omega = omega.sim)
# #Simulate poisson counts
# set.seed(234)
# spatial.df$counts <- rpois(nrow(spatial.df), exp(Beta0+spatial.df$omega))
# spatial.df$mean_counts <- exp(Beta0+spatial.df$omega)
# save(spatial.df, file = "docs/02-intermediate-tmb/lecture_code/spatial_df.RData")

load("docs/02-intermediate-tmb/lecture_code/spatial_df.RData")

sum(spatial.df$mean_counts)
#Plot spatial field
ggplot(spatial.df, aes(x=x, y=y, color=omega)) + geom_point()
#Plot species counts
ggplot(spatial.df, aes(x=x, y=y, color=mean_counts)) + geom_point()

#Sample locations
set.seed(345)
sample.df <- spatial.df[sample(1:nrow(spatial.df), 200),]
Loc <- as.matrix(sample.df[,1:2])

#Create INLA mesh
mesh <- inla.mesh.create(Loc)
mesh <- inla.mesh.2d(Loc, max.edge = 5)

ggplot(sample.df, aes(x=x, y=y)) + geom_point() + gg(mesh) +
  theme_classic()

#Create SPDE object
spde <- inla.spde2.matern(mesh)

# Compile
compile("docs/02-intermediate-tmb/lecture_code/spatial_spde.cpp" )
dyn.load( dynlib("docs/02-intermediate-tmb/lecture_code/spatial_spde") )

## -- Fit TMB model to estimate parameters -- ##
#Set up TMB model
Data <- list(Y_i = sample.df$counts, v_i = mesh$idx$loc-1)
Data$spde <- spde$param.inla[c("M0","M1","M2")]
Params <- list(beta0 = 0, 
               ln_kappa = 0, ln_tau = 0,
               omega_v = rep(0, mesh$n))
obj <- MakeADFun(Data, Params, random = "omega_v", DLL = "spatial_spde")
a <- Sys.time()
opt <- nlminb(obj$par, obj$fn, obj$gr)
b <- Sys.time()

cppad.opt.time <- difftime(b,a)

#Report back results
a <- Sys.time()
sdr <- sdreport(obj)
b <- Sys.time()

cppad.sdr.time <- difftime(b,a)

c(Beta0, log(Kappa), log(Tau))
sdr
report <- obj$report()
plot(sample.df$counts, report$lambda);abline(0,1)

#Map estimated RE field
omega.est <- obj$env$parList()$omega_v
P <- inla.mesh.project(mesh, as.matrix(proj.grid))
A <- inla.spde.make.A(mesh, Loc)
spatial.df$omega_est <- as.vector(P$A %*% omega.est)
spatial.df$mean_est <- as.vector(P$A %*% t(  report$lambda %*% A ))
ggplot(spatial.df, aes(x=x, y=y, color = omega_est)) + geom_point()
ggplot(spatial.df, aes(x=x, y=y, color = omega)) + geom_point()
ggplot(spatial.df, aes(x=x, y=y, color = mean_est)) + geom_point()
ggplot(spatial.df, aes(x=x, y=y, color = mean_counts)) + geom_point()


## Recompile with TMBad and compare speed
dyn.unload(dynlib("docs/02-intermediate-tmb/lecture_code/spatial_spde"))

compile("docs/02-intermediate-tmb/lecture_code/spatial_spde.cpp",
        framework = "TMBad")
dyn.load( dynlib("docs/02-intermediate-tmb/lecture_code/spatial_spde") )
obj <- MakeADFun(Data, Params, random = "omega_v", DLL = "spatial_spde")
a <- Sys.time()
opt <- nlminb(obj$par, obj$fn, obj$gr)
b <- Sys.time()

tmbad.opt.time <- difftime(b,a)

#Report back results
a <- Sys.time()
sdr <- sdreport(obj)
b <- Sys.time()

tmbad.sdr.time <- difftime(b,a)

tmbad.opt.time
cppad.opt.time

tmbad.sdr.time
cppad.sdr.time


## ---- Diagnostics ---- ##

## check laplace approximation 
chk <- checkConsistency(obj)
summary(chk)

## Bias Correction 
# Report out derived mean abundance
summary(sdr, "report")[1,]
est_mean_count <- summary(sdr, "report")[1,]
est_mean_count[1] + c(-1,1)*1.96 * est_mean_count[2]
sum(sample.df$mean_counts)

# Calculate with bias correction
sdr.biascorr <- sdreport(obj, bias.correct = TRUE,
                         bias.correct.control = list(sd = TRUE, split = 1))

summary(sdr.biascorr, "report")[1,]
est_bc_mean_count <- summary(sdr.biascorr, "report")[1,3:4]
est_mean_count[1] + c(-1,1)*1.96 * est_mean_count[2]
est_bc_mean_count[1] + c(-1,1)*1.96 * est_bc_mean_count[2]

## Likelihood Profile
prof <- TMB::tmbprofile(obj, "ln_kappa")
#png(filename = "docs/02-intermediate-tmb/static/ln_kappa_profile.png")
plot(prof);abline(v = log(Kappa))
ci <- confint(prof)
summary(sdr, "fixed")
abline(v = ci[1]); abline(v = ci[2])
abline(v = log(Kappa), col = 'red')
#dev.off()

#Likelihood Profile confidence interval
ci
#Asymptotic confidence interval
summary(sdr, "fixed")[2,1] + c(-1,1) * 1.96 * summary(sdr, "fixed")[2,2]
#pull out omega values at sample locations
omega.est <- as.list(sdr, "Estimate")$omega[mesh$idx$loc]
omega.se <- as.list(sdr, "Std. Error")$omega[mesh$idx$loc]
df <- data.frame(omega.est = omega.est,
                 omega.se = omega.se,
                 conf.lower = omega.est - qnorm(.975)*omega.se,
                 conf.upper = omega.est + qnorm(.975)*omega.se)
#png(filename = "docs/02-intermediate-tmb/static/omega-confint.png")
ggplot(df, aes(x = 1:200, y = omega.est)) +
  geom_point(size = 0.5) +
  geom_ribbon(aes(x = 1:200, ymin = conf.lower, 
                  ymax = conf.upper, fill = "band"), alpha = 0.3) +
  scale_fill_manual("", values = "grey20") +
  ylab("omega") + xlab("") +
  theme_classic() + 
  theme(legend.position = "none")
#dev.off()
