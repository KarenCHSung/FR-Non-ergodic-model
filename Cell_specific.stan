/* *******************************************************************
cell-specific approach
******************************************************************* */
  
  data {
    int<lower=1> N;      // the number of record
    int<lower=1> NCELL;  // the number of Cell
    int<lower=1> NEQ;    // the number of EQ 
    int<lower=1> NSTAT;  // the number of ST
    
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
  real<lower=0> sigma_cA;   // std of cell-specific attenuation
  //mean value
  real <upper=0> mu_cA;       // global (mean) attenuation
  real intercept;
  //coeff
  vector<upper=0>[NCELL] c_Ac; // cell-specific attenuation
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
  mu_cA ~ normal(0,0.01);  // prior for global attenuation parameter
  c_Ac ~ normal(mu_cA,sigma_cA);  // prior for regional attenuation
  
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
