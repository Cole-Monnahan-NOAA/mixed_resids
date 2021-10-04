#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_VECTOR_INDICATOR(keep,y);
  PARAMETER_VECTOR(beta); //intercept
  PARAMETER(ln_sig);
  Type sig = exp(ln_sig);
  Type nll = 0;
  Type cdf;

  vector<Type> mu = X * beta; 
  for(int i=0; i<y.size(); i++){
    nll -= dnorm(y(i), mu(i), sig, true);
    cdf = squeeze(pnorm(y(i),  mu(i), sig ));
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
  }

  SIMULATE{
    y = rnorm(mu, sig);
    REPORT(y);
  }
  
  vector<Type> fpr = mu;
  vector<Type> exp_val = mu;
  REPORT(mu);
  REPORT(exp_val);
  REPORT(fpr);
  REPORT(sig);
  return(nll);
}
