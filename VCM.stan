/*********************************************
This model explicitly estimates the latent (uncorrelated) event terms and station terms
The model does not includes cell-specific attenuation
 ********************************************/

data {
  int N;  // number of records
  int NEQ;  // number of earthquakes
  int NSTAT;  // number of stations
  vector[N] Y; // residual
  //vector[N] mu_rec; 
  int<lower=1,upper=NEQ> eq[N]; // event id (in numerical order from 1 to last)
  int<lower=1,upper=NSTAT> stat[N]; // station id (in numerical order from 1 to last)
  vector[2] X_e[NEQ];  // event coordinate for each record
  vector[2] X_s[NSTAT];  // station coordinate for each record
  vector[NSTAT] lnVS;  //lnVs30

}

transformed data {
  real delta = 1e-9;
}

parameters {
  //real b1;
  real<lower=0> phi_0;  // phi_0 
  real<lower=0> tau_0;  // tau_0 
  real<lower=0> phi_S2S;  // phi_s2s 

  real intercept;
  real <lower=0>rho_eq;
  real <lower=0.1>theta_eq;
  vector<lower=0>[2] rho_stat;
  vector<lower=0.1>[2] theta_stat;

  vector[NEQ] z_eq;
  vector[NSTAT] z_stat;
  vector[NSTAT] z_vs;

  vector[NEQ] eqterm;
  vector[NSTAT] statterm;

}

model {
  vector[NEQ] f_eq;
  vector[NSTAT] f_stat;
  vector[NSTAT] f_vs;
  vector[N] mu_rec2;

  phi_0 ~ lognormal(-1.30,0.3);
  tau_0 ~ lognormal(-1,0.3);
  phi_S2S ~ lognormal(-0.8,0.3);
  
  intercept ~ normal(0,0.1);
  rho_eq ~ inv_gamma(2,0.5);
  rho_stat[1] ~ inv_gamma(2,0.5);
  rho_stat[2] ~ inv_gamma(2,0.5);
  theta_eq ~ exponential(20);
  theta_stat[1] ~ exponential(20);
  theta_stat[2] ~ exponential(20);

  z_eq ~ std_normal();
  z_stat ~ std_normal();
  z_vs ~ std_normal();

  eqterm ~ normal(0,tau_0);
  statterm ~ normal(0,phi_S2S);


  // latent variable event contributions to GP
  {
    matrix[NEQ,NEQ] cov_eq;
    matrix[NEQ,NEQ] L_eq;

    for(i in 1:NEQ) {
      for(j in i:NEQ) {
        real d_e;
        real c_eq;
  
        d_e = distance(X_e[i],X_e[j]);
  
        c_eq = (theta_eq^2 * exp(-d_e/rho_eq));
  
        cov_eq[i,j] = c_eq;
        cov_eq[j,i] = cov_eq[i,j];
      }
      cov_eq[i,i] = cov_eq[i,i] + delta;
    }

    L_eq = cholesky_decompose(cov_eq);
    f_eq = L_eq * z_eq;
  }

  // latent variable station contributions to GP
  { 
    
    matrix[NSTAT,NSTAT] cov_stat;
    matrix[NSTAT,NSTAT] L_stat;
    matrix[NSTAT,NSTAT] cov_vs;
    matrix[NSTAT,NSTAT] L_vs;

    for(i in 1:NSTAT) {
      for(j in i:NSTAT) {
        real d_s;
        real c_vs;
        real c_stat;
  
        d_s = distance(X_s[i],X_s[j]);
        
        c_stat = (theta_stat[1]^2 * exp(-d_s/rho_stat[1]));
        c_vs = lnVS[i] * lnVS[j] * (theta_stat[2]^2 * exp(-d_s/rho_stat[2]));
  
        cov_stat[i,j] = c_stat;
        cov_stat[j,i] = cov_stat[i,j];
        
        cov_vs[i,j] = c_vs;
        cov_vs[j,i] = cov_vs[i,j];
        
      }
      cov_stat[i,i] = cov_stat[i,i] + delta;
      cov_vs[i,i] = cov_vs[i,i] + delta;
    }

    L_stat = cholesky_decompose(cov_stat);
    f_stat = L_stat * z_stat;
    
    L_vs = cholesky_decompose(cov_vs);
    f_vs = L_vs * z_vs;
  }

  mu_rec2 = intercept + f_eq[eq] + f_stat[stat] + f_vs[stat] + eqterm[eq] + statterm[stat];

  Y ~ normal(mu_rec2,phi_0);

}
