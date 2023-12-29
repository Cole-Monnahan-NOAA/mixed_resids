#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_VECTOR_INDICATOR(keep,y);
  PARAMETER_VECTOR(beta); //intercept
  PARAMETER(ln_sig_y);
  Type sig_y = exp(ln_sig_y);
  Type nll = 0;
  Type cdf;

  vector<Type> mu = X * beta; 
  for(int i=0; i<y.size(); i++){
    nll -= dnorm(y(i), mu(i), sig_y, true);
    cdf = squeeze(pnorm(y(i),  mu(i), sig_y ));
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
  }

  SIMULATE{
    y = rnorm(mu, sig_y);
    REPORT(y);
  }
  
  vector<Type> fpr = mu;
  vector<Type> exp_val = mu;
  REPORT(mu);
  REPORT(exp_val);
  REPORT(fpr);
  REPORT(sig_y);
  return(nll);
}
