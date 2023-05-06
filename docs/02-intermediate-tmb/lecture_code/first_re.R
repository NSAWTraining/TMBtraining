require(TMB)
compile("docs/02-intermediate-tmb/lecture_code/first_re.cpp")
dyn.load(dynlib("docs/02-intermediate-tmb/lecture_code/first_re"))

## Test data
set.seed(123)
yr <- rep(1900:2010,each=2)
year <- factor(yr)
quarter <- factor(rep(1:4,length.out=length(year)))
period <- factor((yr > mean(yr))+1)
## Random year+quarter effect, fixed period effect:
B <- model.matrix(~year+quarter-1)
X <- model.matrix(~period-1)
B <- as(B,"dgTMatrix")
X <- as(X,"dgTMatrix")
u <- rnorm(ncol(B)) ## logsdu=0
beta <- rnorm(ncol(X))*100
eps <- rnorm(nrow(B),sd=1) ## logsd0=0
y <- as.numeric( X %*% beta + B %*% u + eps )

#visualize data
plot(X%*%beta, y)
plot(B%*%u, y)

## Fit model
obj <- MakeADFun(data=list(y=y, B=B, X=X),
                 parameters=list(u=u*0, beta=beta*0, lnSDu=1, lnSDy=1),
                 random="u", 
                 DLL="first_re"
                 )
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
beta

report <- obj$report()
report$nll1 #random effect likelihood
report$nll2 #data likelihood 
report$nll #joint likelihood

#marginal likelihood (via Laplace approximation)
obj$fn(opt$par)


sdr <- sdreport(obj, getJointPrecision  = TRUE)
#Calculate Laplace by hand
#Based on the following property of log(det(C))
# C = LL^T, L = lower triangular cholesky
# log(det(C)) = 2trace(log(L))
H <- solve(sdr$cov.fixed)
L <- chol(H)
2*sum(log(diag(L)))

-4/2 * log(2*pi)  + 
  0.5 *log(det(H)) +
  report$nll
H <-  obj$env$spHess(obj$env$last.par.best, random = TRUE)
L <- t(Cholesky(H))
L <- Cholesky(H, perm = TRUE, LDL = FALSE, super = TRUE)
obj$env$h(theta = obj$env$last.par.best, order = 0, hessian = H, L = t(chol(H)))

logdetH <- 2 * determinant(L)$mod
-length(u)/2 * log(2*pi) + 0.5 * logdetH + report$nll

obj$env$f(obj$env$last.par.best, order = 0) + 0.5 * logdetH - length(u)/2 * 
  log(2 * pi)

summary(sdr, "fixed")
summary(sdr, "random")
sdr$pdHess


Hr <- obj$env$spHess(random=TRUE)
Ha <- obj$env$spHess()
Hf1 <- solve(sdr$cov.fixed)
Hf2 <- diag(1/summary(sdr,"fixed")[,2]^2)
2*sum(log(diag(chol(Hr))))
2*sum(log(diag(chol(Hf1))))
2*sum(log(diag(chol(Hf2))))
2*sum(log(diag(chol(Ha))))
log(chol(Hr))
determinant(Hr)$mod

logdetH <- 2 * determinant(Hr)$mod
ans <- f(obj$env$last.par.best, order = 0) + 0.5 * logdetH - length(random)/2 * 
  log(2 * pi)


sdr$cov.fixed
#output sparse Hessian
obj$env$spHess()
