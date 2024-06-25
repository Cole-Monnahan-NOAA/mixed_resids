#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_VECTOR_INDICATOR(keep,y);
  PARAMETER(beta); //intercept
  Type nll = 0;
  Type cdf;

  for(int i=0; i<y.size(); i++){
    nll -= dpois(y(i), exp(beta), true);
    cdf = squeeze(ppois(y(i), exp(beta) ));
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
  }

  SIMULATE{
    y = rpois(exp(beta));
    REPORT(y);
  }
  
  vector<Type> fpr(y.size());
  fpr.fill(beta);
  vector<Type> exp_val(y.size());
  exp_val.fill(beta);
  REPORT(exp_val);
  REPORT(fpr);
  return(nll);
}
