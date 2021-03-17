// Estimate and validate a random walk model with and without drift
//
// Compare Thygesen et al (submitted, 2016): Validation of state space models
// fitted as mixed effects models
//
// Uffe Høgsbro Thygesen and Kasper Kristensen, 2016

// Modifed by Cole to add simulation component for comparing OSA
// vs DHARMa resids. CCM 11/2020

#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);                  // Observations
  DATA_INTEGER(simRE);
  DATA_VECTOR_INDICATOR(keep, y);  // For one-step predictions

  DATA_SCALAR(huge);
  PARAMETER_VECTOR(x);
  PARAMETER(mu);
  PARAMETER(logsigma);
  PARAMETER(logtau);

  Type sigma=exp(logsigma);	// observation sd
  Type tau=exp(logtau);		// process sd
  // Initial condition
  //Type nll = -dnorm(x(0), Type(0), huge, true);
  vector<Type> ypred(x.size());
  ypred.setZero();

  Type nll, nll2;
  nll=0; nll2=0;
  
  // Increments
  if(simRE==1){
    // Simulate new RE if using unconditional simulation resids
    for (int i = 0; i < x.size(); i++){
      SIMULATE{
	x(i)=rnorm(Type(0), tau);
      }
    }
  }
  ypred(0)=x(0); // expected value
  for (int i = 1; i < x.size(); ++i){
    ypred(i)=ypred(i-1)+x(i)+mu;
    // prob of random effects
    nll -= dnorm(x(i), Type(0.0), tau, true); //dnorm(x(i), x(i - 1) + mu, sigma, true);
  } 

  // Data likelihood
  Type cdf;
  for (int i = 0; i < y.size(); ++i){
    nll2 -= keep(i) * dnorm(y(i), ypred(i), sigma, true);
    cdf = squeeze( pnorm(y(i), ypred(i), sigma) );
    nll2 -= keep.cdf_lower(i) * log( cdf );
    nll2 -= keep.cdf_upper(i) * log( 1.0 - cdf );
    SIMULATE {
      y(i) = rnorm(ypred(i), sigma); // conditional
    }
  }

  SIMULATE{
    if(simRE==1){
      REPORT(x);
    }
    REPORT(y);
  }

  REPORT(x);
  REPORT(ypred);
  REPORT(tau);
  REPORT(sigma);
  REPORT(nll);
  REPORT(nll2);
  return nll+nll2;
}
