#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  
  DATA_VECTOR(Codend);
  DATA_VECTOR(Cover);
  DATA_VECTOR(Length);
  
  PARAMETER(beta0);
  PARAMETER(beta1);
  
  Type nll;
  vector<Type> N = Codend + Cover;
  vector<Type> eta  = beta0 + beta1 * Length;
  vector<Type> p = 1/(1+exp(-eta));

  nll = -sum(dbinom(Codend, N, p, true));
  
  Type L50 = -beta0/beta1;
  REPORT(L50);
  ADREPORT(L50);
  REPORT(p);
  
  return nll;
  
}
