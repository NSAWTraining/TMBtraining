library(TMB)

compile('docs/01-beginner-tmb/lecture_code/linReg.cpp',"-O1 -g",DLLFLAGS="")
dyn.load(dynlib('docs/01-beginner-tmb/lecture_code/linReg'))
Data <- list(y = c(-2.1, 3.3, 4.2),
             X = cbind(c(1,1,1),c(4,5,6))
)

Pars <- list(beta = c(0,0), lnSigma = 0)
obj <- MakeADFun(data = Data, 
                 parameters = Pars,  
                 DLL = "linReg")