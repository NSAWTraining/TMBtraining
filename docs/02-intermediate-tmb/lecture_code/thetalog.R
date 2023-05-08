library(TMB)
#use path from docs when running code
compile("docs/02-intermediate-tmb/lecture_code/thetalog.cpp")
dyn.load(dynlib("docs/02-intermediate-tmb/lecture_code/thetalog"))
#use path from 02-intermediate-tmb when generating Rmarkdown
#compile("lecture_code/thetalog.cpp")
#dyn.load(dynlib("lecture_code/thetalog"))

## Read data
Y <- scan("docs/02-intermediate-tmb/lecture_code/thetalog.dat", skip=3, quiet=TRUE)
#Y <- scan("lecture_code/thetalog.dat", skip=3, quiet=TRUE)
data <- list(Y=Y)

## Parameter initial guess
parameters <- list(
  X = data$Y*0,
  lnr0 = 0,
  lntheta = 0,
  lnK = 6,
  lnsdx = 0,
  lnsdy = 0
)

## Fit model
obj <- MakeADFun(data, parameters, random="X", DLL="thetalog")
opt <- nlminb(obj$par, obj$fn, obj$gr)
report <- obj$report()
sdr <- sdreport(obj)
sdr
summary(sdr, "fixed")
summary(sdr, "random")

X.est <- obj$env$parList()$X

# plot(1:length(Y), Y)
# lines(1:length(Y), X.est, col='red', lwd = 2)

report$nll_re
report$nll_obs
report$nll
opt$objective
obj$fn()

# -log marginal likelihood:
# -nRE/2 * log(2*pi) + 0.5 * log(det(H)) + joint nll

# Get Hessian
H <- obj$env$spHess(obj$env$last.par.best, random = TRUE)
# log(det(C)) = 2trace(log(L))
# Use: C = LL^T, L = lower triangular cholesky
L <- chol(H)
logdetH <- 2*sum(log(diag(L)))
#Calculate the Laplace approximation
-length(parameters$X)/2 * log(2*pi) + 0.5 * logdetH + report$nll

opt$objective
obj$fn()

## Run TMB model in Stan using tmbstan
# install.packages(tmbstan)
# library(tmbstan)
# tmbstan(obj)
