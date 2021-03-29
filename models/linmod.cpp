#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_VECTOR_INDICATOR(keep,y);
  PARAMETER(b0); //intercept
  PARAMETER(b1); // slope
  PARAMETER(logsigma);
  Type sigma = exp(logsigma);
  Type nll = 0;
  Type cdf;

  vector<Type> mu(y.size());
  for(int i=0; i<y.size(); i++){
    mu(i) = b0 + b1*x(i);
    nll -= dnorm(y(i), mu(i), sigma, true);
    cdf = squeeze(pnorm(y(i),  mu(i), sigma ));
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
  }

  SIMULATE{
    y = rnorm(mu, sigma);
    REPORT(y);
  }
  
  REPORT(mu);
  REPORT(sigma);
  return(nll);
}
