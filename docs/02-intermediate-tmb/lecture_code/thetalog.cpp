// Theta logistic population model from Pedersen et al 2012, Ecol. Modelling.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data section */
  DATA_VECTOR(Y);
  /* Parameter section */
  PARAMETER_VECTOR(X);
  PARAMETER(lnr0);
  PARAMETER(lntheta);
  PARAMETER(lnK);
  PARAMETER(lnsdx);
  PARAMETER(lnsdy);
  /* Procedure section */
  Type r0 = exp(lnr0);
  Type theta = exp(lntheta);
  Type K = exp(lnK);
  Type sdx = exp(lnsdx);
  Type sdy = exp(lnsdy);
  int timeSteps = Y.size();
 
  Type nll_re = 0;
  Type nll_obs = 0;
  for(int i=1; i<timeSteps; i++){
    Type m = X[i-1] + r0 * (1.0 - pow( exp( X[i-1] ) / K , theta));
    nll_re -= dnorm(X[i], m, sdx, true);
  }
  for(int i=0; i<timeSteps; i++){
    nll_obs -= dnorm(Y[i], X[i], sdy, true);
  }
  
  Type nll = nll_re + nll_obs;
  REPORT(nll_re);
  REPORT(nll_obs);
  REPORT(nll);
  
  return nll;
}
