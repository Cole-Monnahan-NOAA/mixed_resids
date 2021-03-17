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
  DATA_VECTOR_INDICATOR(keep, y);  // For one-step predictions

  DATA_SCALAR(huge);
  PARAMETER_VECTOR(x);
  PARAMETER(mu);
  PARAMETER(logsigma);
  PARAMETER(logs);

  // Initial condition
  Type nll = -dnorm(x(0), Type(0), huge, true);
  // This is for the unconditional DHARMa residuals
  vector<Type> x2(x.size());
  vector<Type> y2(x.size());
  vector<Type> xpred(x.size());
  xpred.setZero();

  // Increments
  SIMULATE{
    x2(0)=rnorm(Type(0) ,exp(logsigma));
  }
  xpred(0)=x(0);
  for (int i = 1; i < x.size(); ++i){
    xpred(i)=xpred(i-1)+x(i)+mu;
    nll -= dnorm(x(i), x(i - 1) + mu, exp(logsigma), true);
    SIMULATE {
      x2(i) = rnorm(x2(i-1)+mu, exp(logsigma));
    }
  }
  // Observations
  Type cdf;
  for (int i = 0; i < y.size(); ++i){
    nll -= keep(i) * dnorm(y(i), x(i), exp(logs), true);
    cdf = squeeze( pnorm(y(i), x(i), exp(logs)) );
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );

    SIMULATE {
      y(i) = rnorm(x(i), exp(logs)); // conditional
      y2(i)=rnorm(x2(i), exp(logs)); // unconditional
    }
  }
  REPORT(x);
  REPORT(y);
  REPORT(x2);
  REPORT(y2);
  REPORT(xpred);
 
  return nll;
}
