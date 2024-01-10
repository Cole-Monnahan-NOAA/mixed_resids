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
  PARAMETER(beta);
  PARAMETER(mu);
  PARAMETER(ln_sig_y);
  PARAMETER_VECTOR(ln_sig_u); 
  PARAMETER_VECTOR(u);

  
  //only fit random effect if ln_sig_u isn't null
  bool u_flag = (ln_sig_u.size()>0);
  Type sig_y = exp(ln_sig_y);	// observation sd
  vector<Type> sig_u(ln_sig_u.size());
  Type nll = 0;
  vector<Type> eta(u.size());
  eta.setZero();
  
  if(u_flag){
    sig_u = exp(ln_sig_u);
    //Partition into deterministic and stochastic components:
    //  y_t = y_0 + at + sum^n_t=1(e_t)
    //Simulate in R using:
    //  z <- cumsum(rnorm(n=N, mean=0, sd=sqrt(variance)))
    //  t <- 1:N
    //  eta <- beta + t*drift + z
    //Initial condition
    nll -= dnorm(u(0), Type(0), Type(1000), true);
    if(sim_re == 1){
      SIMULATE{
        u(0) = rnorm(Type(0), sig_u(0)); 
      }
    }
    for(int i = 1; i < u.size(); i++){
      nll -= dnorm(u(i), u(i-1), sig_u(0), true);
      if(sim_re == 1){
        SIMULATE{
          u(i) = rnorm(u(i-1), sig_u(0));
        }
      }
    }
  }

  // Observations
  Type cdf = 0;
  for (int i = 0; i < y.size(); ++i){
    eta(i) = beta + mu * (i+1) + u(i);
    if(mod == 0){
      nll -= keep(i) * dnorm(y(i), eta(i), sig_y, true);
      cdf = squeeze( pnorm(y(i), eta(i), sig_y) );
      SIMULATE {
        y(i) = rnorm(eta(i), sig_y); // conditional
      }
    }
    if(mod == 1){
      nll -= keep(i) * dlnorm(y(i), eta(i), sig_y, true);
      cdf = squeeze( pnorm(log(y(i)), eta(i), sig_y) );
      SIMULATE {
        y(i) = exp(rnorm(eta(i), sig_y)); // conditional
      }
    }
    if(mod == 2){ 
      nll -= keep(i) * dgamma( y(i), 1/pow(sig_y,2), exp(eta(i))*pow(sig_y,2), true);
      cdf = squeeze( pgamma(y(i), 1/pow(sig_y,2), exp(eta(i))*pow(sig_y,2)) );
      SIMULATE{
        y(i) = rgamma( 1/pow(sig_y,2), exp(eta(i))*pow(sig_y,2) );
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
  for(int i=0; i<y.size(); i++) fpr(i) =  beta + mu * (i+1);
  vector<Type> exp_val(y.size());
  if(mod == 2){
    exp_val = exp(eta);
    fpr = exp(fpr);
  } else {
    exp_val = eta;
  } 
  
  REPORT(u);
  REPORT(exp_val);
  REPORT(fpr);
  REPORT(eta);
  REPORT(sig_u);
  REPORT(sig_y);
  REPORT(nll);
  
  return nll;
}
