#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);
  DATA_IVECTOR(group);
  DATA_INTEGER(sim_re);
  PARAMETER(b0); //intercept
  PARAMETER(ln_sig_u); //natural log of between subject sd 
  PARAMETER(ln_sig_y); 
  PARAMETER_VECTOR(u); //random effects

  DATA_VECTOR_INDICATOR(keep,y);

  Type sig_u = exp(ln_sig_u);
  Type sig_y = exp(ln_sig_y);
  Type nll = 0;

  Type cdf;

  vector<Type> mu(y.size());
  for(int j=0; j<u.size(); j++){
    nll -= dnorm(u(j), Type(0), sig_u, true);
    if(sim_re == 1){
       SIMULATE{
       u(j) = rnorm(Type(0), sig_u);
     }
    }
    
    for(int i=0; i<y.size(); i++){
      mu(i) = b0 + u(group(i));
      nll -= keep(i) * dnorm(y(i), mu(i), sig_y, true);
      cdf = squeeze(pnorm(y(i),  mu(i), sig_y ));
      nll -= keep.cdf_lower(i) * log( cdf );
      nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
      SIMULATE{
        y(i) = rnorm(mu(i), sig_y);
      }
    }
  }
  SIMULATE{
    REPORT(y);
    if(sim_re == 1){
      REPORT(u);
    }
  }
  REPORT(mu);
  REPORT(u);
  REPORT(sig_u);
  REPORT(sig_y);

  return(nll);
}
