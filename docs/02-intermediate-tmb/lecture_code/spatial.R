# package dependencies: mvtnorm
library(ggplot2)
library(viridis)

# Matern correlation function 
cMatern <- function(H, Nu, R) {
  ifelse(H > 0, besselK(H/R, Nu) * (H/R)^Nu, 1) / gamma(Nu) * 2^(1-Nu)
}

## Simulate Data========================================

Dim = c("n_x"=20, "n_y"=20)
loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])
# distance matrix
Dmat <- as.matrix(dist(loc_xy), diag = TRUE, upper = TRUE)

beta0 <- 1
Phi <- 1
Power <- 1.4
Range <- 2 
sig2 <- 2

## Simulate spatial process
Sigma <- sig2 * cMatern(Dmat,1,Range)
plot(as.vector(Dmat), as.vector(Sigma))

set.seed(123)
omega <- t(mvtnorm::rmvnorm(1, rep(0,nrow(Sigma)), 
                        sigma = Sigma, method = 'chol'))

df <- data.frame(x = loc_xy$x, y = loc_xy$y, omega = as.vector(omega))
ggplot(df, aes(x=x, y=y, fill = omega)) + geom_tile() +
  scale_fill_gradient2() + theme_classic()

## Simulate Tweedie weights
mu <- as.vector(exp(beta0 + omega))
y <- tweedie::rtweedie(nrow(loc_xy), mu = mu, phi = Phi, power = Power)
hist(y, breaks = 100)
df$weights <- y
ggplot(df, aes(x=x, y=y, fill = weights)) + geom_tile() +
  scale_fill_viridis()

total_weight =  sum( mu )

## Fit Model ==================================
library(TMB)
#compile('2022_TMB_Session_II/src/spatial.cpp', framework = "TMBad")
dyn.load(dynlib('2022_TMB_Session_II/src/spatial'))

Data <- list(
  y = y, 
  D = Dmat
)

Pars <- list(
  beta0 = 0,
  ln_phi = 0,
  logit_power = 0,
  ln_range = 0,
  ln_sigma2 = 0,
  omega = rep(0, length(y))
)
obj <- TMB::MakeADFun( 
  data=Data, 
  parameters=Pars, 
  random="omega", DLL="spatial" ,
  hessian = TRUE)
# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)
report <- obj$report()


## check laplace approximation ==================================
chk <- checkConsistency(obj)
summary(chk)

## Bias Correction ==============================================

# Calculate SEs
sdr <- sdreport(obj)
summary(sdr,"fixed")
cbind(summary(sdr, "report"), total_weight)

# Calculate with bias correction
sdr <- sdreport(obj, bias.correct = TRUE,
                bias.correct.control = list(sd = TRUE))
cbind(summary(sdr, "report"), total_weight)

#Plot Confidence intervals
library(ggplot2)
est <- summary(sdr, "report")[c(1,3)]
conf.lower <- est - 1.96*summary(sdr, "report")[c(2,4)]
conf.upper <- est + 1.96*summary(sdr, "report")[c(2,4)]

df <- data.frame(est = est, 
                 lower = conf.lower, 
                 upper = conf.upper, 
                 type = c("w/o bias correction", "bias correction"),
                 true = c("true", "true")
                 )
ggplot(df, aes(x = est, y = type)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower, xmax = upper, y = type)) +
  geom_vline(aes(xintercept = total_weight, color = true), show.legend = TRUE)+ 
  scale_color_manual(name = "", values = c(true = "red")) +
  theme_classic()


## Likelihood Profile

prof <- TMB::tmbprofile(obj, "ln_sigma2")
plot(prof)
ci <- confint(prof)
summary(sdr, "fixed")
abline(v = ci[1]); abline(v = ci[2])
abline(v = log(sig2), col = 'red')

#Likelihood Profile confidence interval
ci
#Asymptotic confidence interval
summary(sdr, "fixed")[5,1] + c(-1,1) * 1.96 * summary(sdr, "fixed")[5,2]

df <- data.frame(omega.est = as.list(sdr, "Estimate")$omega,
                 omega.se = as.list(sdr, "Std. Error")$omega,
                 conf.lower = omega.est - qnorm(.975)*omega.se,
                 conf.upper = omega.est + qnorm(.975)*omega.se)

png(filename = "2022_TMB_Session_II/docs/static/omega-confint.png")
ggplot(df, aes(x = 1:400, y = omega.est)) +
  geom_point() +
  geom_ribbon(aes(x = 1:400, ymin = conf.lower, 
                  ymax = conf.upper, fill = "band"), alpha = 0.3) +
  scale_fill_manual("", values = "grey20") +
  ylab("omega") + xlab("") +
  theme_classic() + 
  theme(legend.position = "none")
dev.off()
