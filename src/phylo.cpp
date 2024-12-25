#include <TMB.hpp>


//see pnbinom_mu at: https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/nmath/pnbinom.c
template<class Type>
Type pnbinom(Type x, Type mu, Type size) {
  Type pr = size/(size + mu);
  return pbeta(pr, size, x+1);
}

template <class Type>
  Type objective_function<Type>::operator()()
{
  using namespace density;
  using namespace Eigen;
    
  DATA_VECTOR(y);     
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(Q);            
  DATA_INTEGER(mod); //0: normal; 1: negative binomial; 2: poisson
  DATA_INTEGER(sim_re);
  DATA_VECTOR_INDICATOR(keep, y);  // For one-step predictions
  
  PARAMETER_VECTOR(beta);
  PARAMETER(theta);
  PARAMETER_VECTOR(ln_sig_u); 
  PARAMETER_VECTOR(u);
  
  // only fit random effect if ln_sig_u isn't null
  bool u_flag = (ln_sig_u.size()>0);
  int n = y.size();
  
  vector<Type> sig_u(ln_sig_u.size());
  Type nll = 0;
  
  if(u_flag){
    sig_u(0) = exp(ln_sig_u(0)); 
    nll += SCALE( GMRF(Q), sig_u(0) )(u);   
    if(sim_re == 1){
      SIMULATE{
        GMRF(Q).simulate(u);
        u = u * sig_u(0);
        REPORT(u);
      }
    }
    REPORT(sig_u);
  }
  
  vector<Type> eta = X*beta;
  Type cdf;
  vector<Type> mu(eta.size());
  Type sig_y = 0;
  Type size = 0;
  Type s1 = 0;
  Type s2 = 0;
  if(mod == 0){
    sig_y = exp(theta);
  }
  if(mod == 1){
    size = exp(theta);
  }
  //Data Likelihood
  for(int i=0; i<n; i++){    
    eta(i) += u(i);
    
    if(mod == 0){ // Gaussian
      mu(i) = eta(i);
      nll -= keep(i) * dnorm(y(i), mu(i), sig_y, true );
      cdf = squeeze( pnorm(y(i), mu(i), sig_y) );
      nll -= keep.cdf_lower(i) * log( cdf );
      nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
      SIMULATE{
        y(i) = rnorm(mu(i), sig_y);
      }
    }
    
    if(mod == 1) { // Neg-Binom
      //variance = mu(1 + mu/size)
      //use dnbinom_robust(y, log(mu), log(var-mu))
      //log(size) = 2 x log(mu) - log(var-mu)
      //log(var-mu) = 2 x log(mu) - log(size)
      mu(i) = exp(eta(i));
      s1 = eta(i); // log(mu_i)
      s2  = 2. * s1 - theta; // log(var - mu) 
      nll -= keep(i) * dnbinom_robust(y(i), s1, s2, true);
      cdf = squeeze( pnbinom(y(i), mu(i), size) );
      nll -= keep.cdf_lower(i) * log( cdf );
      nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
      SIMULATE { // from glmmTMB: uses rnbinom2(mu, var)
        s1 = mu(i);
        s2 = mu(i) * (Type(1) + mu(i) / size);
        y(i) = rnbinom2(s1, s2);
      }
    }
    
    if(mod == 2) { //Poisson
      mu(i) = exp(eta(i));
      nll -= keep(i) * dpois(y(i), mu(i), true);
      cdf = squeeze( ppois(y(i), mu(i)) );
      nll -= keep.cdf_lower(i) * log( cdf );
      nll -= keep.cdf_upper(i) * log( 1.0 - cdf );
      SIMULATE{
        y(i) = rpois(mu(i));
      }
    }
  }
    
  SIMULATE{
    REPORT(y);
  }
  
  vector<Type> fpr = X * beta;
  vector<Type> exp_val = mu;
  
  REPORT(eta);
  REPORT(mu);
  REPORT(fpr);
  REPORT(exp_val);
  REPORT(sig_y);
  REPORT(nll);
  
  return nll;
  
}
