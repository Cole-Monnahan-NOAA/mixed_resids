#include <TMB.hpp>

enum valid_family{
  gaussian_family = 000,
  Gamma_family = 100,
  Poisson_family = 200,
  lognormal_family = 300,
  Tweedie_family = 400,
  Delta_Gamma_family = 500,
  NB_family = 600
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
  PARAMETER_VECTOR(ln_sig_y); //tweedie: this is phi parameter, var = mu^pow/phi
  PARAMETER(theta) //tweedie: logit power; nb: log size; delta: logit prob zero
  PARAMETER_VECTOR(ln_sig_u); //natural log of between subject sd  
  PARAMETER_VECTOR(u); //group-level random effects

  DATA_VECTOR_INDICATOR(keep,y);

  Type sig_y = 0;
  Type sig_u = 0;
  Type power = 0;
  Type pz = 0;
  Type size = 0;
  if(ln_sig_y.size()>0) sig_y = exp(ln_sig_y(0));

  if( family == 400 ){ //Tweedie
    power = invlogit(theta) + Type(1); //tweedie power parameter
  }
  if( family == 500 ){ //Delta
    pz = invlogit(theta); //prob zero in delta model
  }
  if(family == 600){ //NB
    size = exp(theta);
  }
  
  Type nll = 0;
  //only fit random effect if ln_sig_u isn't null
  bool u_flag = (ln_sig_u.size()>0);

  Type cdf;

  // Probability of group-level random effects
  if(u_flag){
    sig_u = exp(ln_sig_u(0));
    for(int j=0; j<u.size(); j++){
      nll -= dnorm(u(j), Type(0), sig_u, true);
      if(sim_re == 1){
        SIMULATE{
	        u(j) = rnorm(Type(0), sig_u);
        }
      }
    }
  }
  
  
  vector<Type> eta = X*beta;
  vector<Type> fpr(y.size());
  vector<Type> mu(y.size());
  Type s1, s2;
  // Likelihood
  for(int i=0; i<y.size(); i++){
    fpr(i) = eta(obs(i));
    fpr(i) = inverse_linkfun(fpr(i), link);
    mu(i) = eta(obs(i)) + u(group(i));
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
      case Tweedie_family:
        nll -= keep(i) * dtweedie(y(i), mu(i), sig_y, power, true);
        //cdf method not possible as no qtweedie in TMB
        SIMULATE{
          y(i) = rtweedie(mu(i), sig_y, power);
        }
        break;
      case NB_family: //variance = mu(1 + mu/size)
        //use dnbinom_robust(y, log(mu), log(var-mu))
        //log(size) = 2 x log(mu) - log(var-mu)
        //log(var-mu) = 2 x log(mu) - log(size)
        s1 = log(mu(i)); // log(mu_i)
        s2  = 2. * s1 - theta; // log(var - mu) //s2 is m (size); 
        nll -= keep(i) * dnbinom_robust(y(i), s1, s2, true);
        //cdf method not possible as no qnbinom in TMB
        SIMULATE { // from glmmTMB: uses rnbinom2(mu, var)
          s1 = mu(i);
          s2 = mu(i) * (Type(1) + mu(i) / exp(theta));
          y(i) = rnbinom2(s1, s2);
        }
        break;
      case Delta_Gamma_family:
        if(y(i) == 0){
          //log(pz)
          nll -= keep(i) * dbinom(Type(0), Type(1), 1-pz, true);
        }
        if(y(i) > 0){
          //log(1-pz) + log(dgamma())
          nll -= keep(i) * ( dbinom(Type(1), Type(1), 1-pz, true) +
            //shape = 1/CV^2; scale = mean*CV^2
            dgamma(y(i), 1/pow(sig_y,2), mu(i)*pow(sig_y,2), true) );
        }
        SIMULATE{
          y(i) = rbinom(Type(1), 1-pz) * 
            rgamma(1/pow(sig_y,2), mu(i)*pow(sig_y,2));
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
  REPORT(sig_y);
  REPORT(sig_u);
  REPORT(power);
  REPORT(pz);
  REPORT(size);

  return(nll);
}
