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
  // Increments
  SIMULATE{x(0)=rnorm(Type(0) ,exp(logsigma));}
  for (int i = 1; i < x.size(); ++i){
    nll -= dnorm(x(i), x(i - 1) + mu, exp(logsigma), true);
    SIMULATE {
      x(i) = rnorm(x(i-1)+mu, exp(logsigma));
    }
  }
  // Observations
  for (int i = 0; i < y.size(); ++i){
    nll -= keep(i) * dnorm(y(i), x(i), exp(logs), true);
    SIMULATE {
      y(i) = rnorm(x(i), exp(logs));
    }
  }
  REPORT(x);
  REPORT(y);
    
  return nll;
}
