#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(obs);
  DATA_IVECTOR(group);
  PARAMETER_VECTOR(beta); //intercepts
  PARAMETER(logsigma);
  vector<Type> mu(obs.size());
  for(int i=0; i<obs.size(); i++)  mu(i)=beta(group(i)-1);
  Type nll = -dnorm(obs, mu, exp(logsigma), true).sum();
  REPORT(mu);REPORT(beta);
  return(nll);
}
