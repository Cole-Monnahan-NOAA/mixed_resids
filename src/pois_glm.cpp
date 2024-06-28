#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_VECTOR_INDICATOR(keep,y);
  PARAMETER(beta); //intercept
  Type nll = 0;
  Type cdf;
  Type lambda = exp(beta);

  for(int i=0; i<y.size(); i++){
    nll -= dpois(y(i), lambda, true);
    cdf = squeeze(ppois(y(i), lambda ));
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
  }

  SIMULATE{
    y = rpois(lambda);
    REPORT(y);
  }
  
  vector<Type> fpr(y.size());
  fpr.fill(lambda);
  vector<Type> exp_val(y.size());
  exp_val.fill(lambda);
  REPORT(exp_val);
  REPORT(fpr);
  return(nll);
}
