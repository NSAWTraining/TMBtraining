// Simple linear regression
#include <TMB.hpp>
template <class Type>
  
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta);
  PARAMETER(lnSigma);
  
  Type sigma = exp(lnSigma);
  
  Type nll = 0;
  int n = y.size();
  matrix<Type> mu = X * beta;
  
  for(int i=0; i<n; i++){
    nll -= dnorm(y(i), mu(i), sigma, true);
  }
  
  Type sig2 = pow(sigma,2);
  REPORT(sig2);
  ADREPORT(sig2);
  
  return nll;
}
