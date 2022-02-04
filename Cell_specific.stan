/* ****************************************************************************
cell-specific approach (Kuehn et al., 2019), the orgignal code is edit by Nico.

Modifications:
(1) For the French dataset, the cell size of RC is  0.2 times 0.2 degrees.
*******************************************************************************/
  
  data {
    int<lower=1> N;      
    int<lower=1> NCELL;  
    int<lower=1> NEQ;    
    int<lower=1> NSTAT;  
    
    vector[N] Y;
    vector[N] mu;
    vector[NCELL] RC[N];
    int<lower=1,upper=NEQ> eq[N];
    int<lower=1,upper=NSTAT> stat[N];
    
  }

parameters {
  //sigma
  real<lower=0> sigma_rec;
  real<lower=0> sigma_eq;
  real<lower=0> sigma_stat;
  real<lower=0> sigma_cA;  
  //mean value
  real <upper=0> mu_cA; 
  real intercept;
  //coeff
  vector<upper=0>[NCELL] c_Ac; 
  //residual
  vector[NEQ] eqterm;
  vector[NSTAT] statterm;
}

model {
  vector[N] mu_rec;
  
  //sigma distribution
  sigma_rec ~ cauchy(0,0.5);
  sigma_cA ~ cauchy(0,0.01);
  // anelastic attenuation
  mu_cA ~ normal(0,0.01);  
  c_Ac ~ normal(mu_cA,sigma_cA);  
  
  //residuals 
  sigma_eq ~ cauchy(0.,0.5);
  sigma_stat ~ cauchy(0.,0.5);
  eqterm ~ normal(0,sigma_eq);
  statterm ~ normal(0,sigma_stat);
  intercept ~ normal(0,0.1);
  
  for(i in 1:N) {
    
    mu_rec[i] = intercept + mu[i] + dot_product(c_Ac,RC[i]) + eqterm[eq[i]] + statterm[stat[i]];
  }
  
  Y ~ normal(mu_rec,sigma_rec);
  
}
