#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  DATA_VECTOR(Y_i); //counts
  DATA_IVECTOR(v_i); //maps mesh vertices to sample locations

  DATA_STRUCT(spde, spde_t); //INLA M0,M1,M2

  PARAMETER(beta0); //intercept
  PARAMETER(ln_alpha); //gamma shape parameter
  PARAMETER(ln_kappa); 
  PARAMETER(ln_tau); 

  //Random Effects vector
  PARAMETER_VECTOR(omega_v);

  Type alpha = exp(ln_alpha);
  Type kappa = exp(ln_kappa);
  Type tau = exp(ln_tau);

  int n_i = Y_i.size();
  int n_v = omega_v.size();

  //Initiate negative log likelihood function
  vector<Type>nll_comp(2);
  nll_comp.setZero();

  //Likelihood of spatial random effects
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  nll_comp(0) = SCALE(GMRF(Q), 1/tau)(omega_v);
  SIMULATE{
    vector<Type>omega(n_v);
    GMRF(Q).simulate(omega);
    omega_v = omega/tau;
  }

  //Likelihood of observations
  vector<Type>ln_mean(n_i);
  for(int i=0; i<n_i; i++){
    ln_mean(i) = beta0 + omega_v(v_i(i));
    nll_comp(1) -= dgamma(Y_i(i), alpha, exp(ln_mean(i))/alpha, true);
    SIMULATE{
      Y_i(i) = rgamma(alpha, exp(ln_mean(i))/alpha);
    }
  }
  
  Type nll = nll_comp.sum();
  Type sig2 = 1 / (4*M_PI*kappa*kappa*tau*tau);

  SIMULATE{
    REPORT(omega_v);
    REPORT(Y_i);
  }
  
  REPORT(ln_mean);
  REPORT(omega_v);
  REPORT(sig2);

  return nll;
}
  
