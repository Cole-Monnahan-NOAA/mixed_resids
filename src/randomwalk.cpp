// Estimate and validate a random walk model with and without drift
//
// Compare Thygesen et al (submitted, 2016): Validation of state space models
// fitted as mixed effects models
//
// Uffe Høgsbro Thygesen and Kasper Kristensen, 2016

#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);                  // Observations
  DATA_INTEGER(sim_re);
  DATA_VECTOR_INDICATOR(keep, y);  // For one-step predictions

  //DATA_SCALAR(huge);
  PARAMETER_VECTOR(u);
  PARAMETER(mu);
  PARAMETER(ln_sig);
  PARAMETER(ln_tau);

  Type sig=exp(ln_sig);	// observation sd
  Type tau=exp(ln_tau);		// process sd
  // Initial condition
  Type nll = -dnorm(u(0), Type(0), Type(1000), true);

  // Increments
  vector<Type> ypred(u.size());
  ypred.setZero();
  ypred(0) = u(0);
  for (int i = 1; i < u.size(); ++i){
    ypred(i) = u(i - 1) + mu;
    nll -= dnorm(u(i), ypred(i), tau, true);
    if(sim_re == 1){
      SIMULATE{
        u(i) = rnorm(ypred(i), tau);
      }
    }
  }

  // Observations
  Type cdf;
  for (int i = 0; i < y.size(); ++i){
    nll -= keep(i) * dnorm(y(i), u(i), sig, true);
    cdf = squeeze( pnorm(y(i), u(i), sig) );
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
    SIMULATE {
      y(i) = rnorm(u(i), sig); // conditional
    }
  }
  
  SIMULATE{
    if(sim_re==1){
      REPORT(u);
    }
    REPORT(y);
  }
  
  vector<Type> fpr(y.size());
  for(int i=0; i<y.size(); i++) fpr(i) = mu;
  vector<Type> exp_val = u;
  
  REPORT(u);
  REPORT(exp_val);
  REPORT(fpr);
  REPORT(tau);
  REPORT(sig);
  REPORT(nll);
  
  return nll;
}
