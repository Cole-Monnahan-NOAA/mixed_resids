// Estimate and validate a random walk model with and without drift
//
// Compare Thygesen et al (submitted, 2016): Validation of state space models
// fitted as mixed effects models
//
// Uffe H?gsbro Thygesen and Kasper Kristensen, 2016

// Modifed by Cole to add simulation component for comparing OSA
// vs DHARMa resids. CCM 11/2020

#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);                  // Observations
  DATA_INTEGER(sim_re);
  DATA_VECTOR_INDICATOR(keep, y);  // For one-step predictions

  PARAMETER_VECTOR(u);
  PARAMETER(mu);
  PARAMETER(ln_sig);
  PARAMETER(ln_tau);

  Type sig=exp(ln_sig);	// observation sd
  Type tau=exp(ln_tau);		// process sd
  // Initial condition
  //Type nll = -dnorm(x(0), Type(0), huge, true);
  vector<Type> ypred(u.size());
  ypred.setZero();

  Type nll, nll2;
  nll=0; nll2=0;
  
  // Increments
  if(sim_re==1){
    // Simulate new RE if using unconditional simulation resids
    for (int i = 0; i < u.size(); i++){
      SIMULATE{
	      u(i)=rnorm(Type(0), tau);
      }
    }
  }
  ypred(0)=u(0); // expected value
  for (int i = 1; i < u.size(); ++i){
    ypred(i)=ypred(i-1)+u(i)+mu;
    // prob of random effects
    nll -= dnorm(u(i), Type(0.0), tau, true); //dnorm(u(i), u(i - 1) + mu, sig, true);
  } 

  // Data likelihood
  Type cdf;
  for (int i = 0; i < y.size(); ++i){
    nll2 -= keep(i) * dnorm(y(i), ypred(i), sig, true);
    cdf = squeeze( pnorm(y(i), ypred(i), sig) );
    nll2 -= keep.cdf_lower(i) * log( cdf );
    nll2 -= keep.cdf_upper(i) * log( 1.0 - cdf );
    SIMULATE {
      y(i) = rnorm(ypred(i), sig); // conditional
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
  vector<Type> exp_val = ypred;

  REPORT(u);
  REPORT(exp_val);
  REPORT(fpr);
  REPORT(tau);
  REPORT(sig);
  REPORT(nll);
  REPORT(nll2);
  return nll+nll2;
}
