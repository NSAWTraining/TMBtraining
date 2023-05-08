// Simple linear regression
#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  
  PARAMETER_VECTOR(beta);
  PARAMETER(lnSigma);
  
  Type nll = 0;
  vector<Type> mu = X*beta;
  Type sigma = exp(lnSigma);

  for(int i=0; i<y.size(); i++){
    nll -= dnorm(y(i), mu(i), sigma, true);
  }
  
  REPORT(sigma);
  ADREPORT(sigma);
  
  return nll;
}
