#include <TMB.hpp>

enum valid_family{
  gaussian_family = 000,
  Gamma_family = 100,
  Poisson_family = 200,
  lognormal_family = 300,
  Tweedie_family = 400
};

enum  valid_link{
  log_link = 0,
  logit_link = 1,
  identity_link = 2, 
  probit_link = 3
};

template<class Type>
Type inverse_linkfun(Type eta, int link){
  Type ans;
  switch(link){
  case log_link:
    ans = exp(eta);
    break;
  case identity_link:
    ans = eta;
    break;
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  default:
    error("Link not implemented");
  } // end switch
  return ans; 
}

template<class Type>
Type objective_function<Type>::operator()()
{
  
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_IVECTOR(group);
  DATA_IVECTOR(obs);
  DATA_INTEGER(sim_re);
  DATA_INTEGER( family );
  DATA_INTEGER( link );
  PARAMETER_VECTOR(beta); //intercept
  PARAMETER_VECTOR(ln_sig_y); 
  PARAMETER(ln_sig_u); //natural log of between subject sd 
  PARAMETER_VECTOR(ln_sig_v); //natural log of overdispersion sd 
  PARAMETER_VECTOR(u); //group-level random effects
  PARAMETER_VECTOR(v); //observation-level random effects

  DATA_VECTOR_INDICATOR(keep,y);

  Type sig_u = exp(ln_sig_u);
  Type sig_y;
  if(ln_sig_y.size()>0)  sig_y = exp(ln_sig_y(0));
  Type nll = 0;
  //only fit overdispersion if ln_sig_v isn't null
  bool v_flag = (ln_sig_v.size()>0);
  Type sig_v= 0;

  Type cdf;

  // Probability of group-level random effects
  for(int j=0; j<u.size(); j++){
    nll -= dnorm(u(j), Type(0), sig_u, true);
    if(sim_re == 1){
      SIMULATE{
	      u(j) = rnorm(Type(0), sig_u);
      }
    }
  }
  // Probability of observation-level random effects
  if(v_flag){
    sig_v = exp(ln_sig_v(0));
    for(int j=0; j<v.size(); j++){
      nll -= dnorm(v(j), Type(0), sig_v, true);
      if(sim_re == 1){
        SIMULATE{
          v(j) = rnorm(Type(0), sig_v);
        }
      }
    }
  }
  
  vector<Type> eta = X*beta;
  vector<Type> fpr(y.size());
  vector<Type> mu(y.size());
  // Likelihood
  for(int i=0; i<y.size(); i++){
    fpr(i) = eta(obs(i));
    fpr(i) = inverse_linkfun(fpr(i), link);
    mu(i) = eta(obs(i)) + u(group(i)) + v(obs(i));
    mu(i) = inverse_linkfun(mu(i), link);

    switch(family){
      case gaussian_family:
        nll -= keep(i) * dnorm(y(i), mu(i), sig_y, true);
        cdf = squeeze( pnorm(y(i), mu(i), sig_y) );
        nll -= keep.cdf_lower(i) * log( cdf );
        nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
        SIMULATE{
          y(i) = rnorm(mu(i), sig_y);
        }
        break;
      case Gamma_family:
        //shape = 1/CV^2; scale = mean*CV^2
        nll -= keep(i) * dgamma( y(i), 1/pow(sig_y,2), mu(i)*pow(sig_y,2), true);
        cdf = squeeze( pgamma(y(i), 1/pow(sig_y,2), mu(i)*pow(sig_y,2)) );
        nll -= keep.cdf_lower(i) * log( cdf );
        nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
        SIMULATE{
          y(i) = rgamma( 1/pow(sig_y,2), mu(i)*pow(sig_y,2) );
        }
        break;
      case Poisson_family:
        nll -= keep(i) * dpois(y(i), mu(i), true);
        cdf = squeeze( ppois(y(i), mu(i)) );
        nll -= keep.cdf_lower(i) * log( cdf );
        nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
        SIMULATE{
          y(i) = rpois(mu(i));
        }
        break;
      default:
        error("Family not supported");
    }
    
  }
  
  SIMULATE{
    REPORT(y);
    if(sim_re == 1){
      REPORT(u);
    }
  }
  
  vector<Type> exp_val = mu;

  REPORT(mu);
  REPORT(fpr);
  REPORT(exp_val);
  REPORT(u);
  REPORT(v);
  REPORT(sig_y);
  REPORT(sig_u);
  REPORT(sig_v);

  return(nll);
}
