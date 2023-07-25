require(TMB)
compile("docs/01-beginner-tmb/lecture_code/first_re.cpp")
dyn.load(dynlib("docs/01-beginner-tmb/lecture_code/first_re"))

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
sdr <- sdreport(obj)

#Calculate Laplace by hand
# -log marginal likelihood:
# -nRE/2 * log(2*pi) + 0.5 * log(det(H)) + joint nll

# Get Hessian
H <-  obj$env$spHess(obj$env$last.par.best, random = TRUE)
# log(det(C)) = 2trace(log(L))
# Use: C = LL^T, L = lower triangular cholesky

L <- chol(H)
logdetH <- 2*sum(log(diag(L)))
#Calculate the Laplace approximation
-length(u)/2 * log(2*pi)  + 0.5 * logdetH + report$nll
opt$objective
obj$fn()

summary(sdr, "fixed")
summary(sdr, "random")
sdr$pdHess

as.list(sdr, "Estimate")$u
as.list(sdr, "Std")$u


