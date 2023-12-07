// Estimate and validate a random walk model with and without drift
//
// Compare Thygesen et al (submitted, 2016): Validation of state space models
// fitted as mixed effects models
//
// Uffe Høgsbro Thygesen and Kasper Kristensen, 2016


#include <TMB.hpp>

template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);                  // Observations
  DATA_INTEGER(mod); //0: normal; 1: lognormal
  DATA_INTEGER(sim_re);
  DATA_VECTOR_INDICATOR(keep, y);  // For one-step predictions

  //DATA_SCALAR(huge);
  PARAMETER_VECTOR(u);
  PARAMETER(mu);
  PARAMETER(ln_sig_y);
  PARAMETER_VECTOR(ln_sig_u); 

  
  //only fit random effect if ln_sig_u isn't null
  bool u_flag = (ln_sig_u.size()>0);

  Type sig_y = exp(ln_sig_y);	// observation sd

  if(u_flag){
    vector<Type> sig_u = exp(ln_sig_u);		// process sd
    // Initial condition
    Type nll = -dnorm(u(0), Type(0), Type(1000), true);
    /* Turn off simulation of initial condition
    if(sim_re == 1){
      SIMULATE{
        u(0) = rnorm(Type(0), sig_u); 
      }
    }
    */
    // Increments
    vector<Type> ypred(u.size());
    ypred.setZero();
    ypred(0) = u(0);

    for (int i = 1; i < u.size(); ++i){
      ypred(i) = u(i - 1) + mu;
      nll -= dnorm(u(i), ypred(i), sig_u, true);
      if(sim_re == 1){
        SIMULATE{
          u(i) = rnorm(ypred(i), sig_u);
        }
      }
    }
  }

  // Observations
  Type cdf = 0;
  for (int i = 0; i < y.size(); ++i){
    if(mod == 0){
      nll -= keep(i) * dnorm(y(i), u(i), sig_y, true);
      cdf = squeeze( pnorm(y(i), u(i), sig_y) );
      SIMULATE {
        y(i) = rnorm(u(i), sig); // conditional
      }
    }
    if(mod == 1){
      nll -= keep(i) * dlnorm(y(i), u(i), sig_y, true);
      cdf = squeeze( pnorm(log(y(i)), u(i), sig_y) );
      SIMULATE {
        y(i) = exp(rnorm(u(i), sig)); // conditional
      }
    }
    if(mod == 2){ 
      nll -= keep(i) * dgamma( y(i), 1/pow(sig_y,2), exp(u(i))*pow(sig_y,2), true);
      cdf = squeeze( pgamma(y(i), 1/pow(sig_y,2), exp(u(i))*pow(sig_y,2)) );
      SIMULATE{
        y(i) = rgamma( 1/pow(sig_y,2), exp(u(i))*pow(sig_y,2) );
      }
    }
    nll -= keep.cdf_lower(i) * log( cdf );
    nll -= keep.cdf_upper(i) * log( 1.0 - cdf );

  }
  
  SIMULATE{
    if(u_flag){
      if(sim_re==1){
        REPORT(u);
      }
    }
    REPORT(y);
  }
  
  vector<Type> fpr(y.size());
  for(int i=0; i<y.size(); i++) fpr(i) = mu;
  vector<Type> exp_val = u;
  
  REPORT(u);
  REPORT(exp_val);
  REPORT(fpr);
  if(u_flag){
    REPORT(ypred);
  }
  REPORT(sig_u);
  REPORT(sig_y);
  REPORT(nll);
  
  return nll;
}
