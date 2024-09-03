// Stan program for regression of the pinene problem
functions{
vector ode_pinene(real t, 
               vector y,
               real k1,
               real k2,
               real k3,
               real k4,
               real k5){
  
  vector[5] dydt; // number of reactions
  
  dydt[1] = - k1*y[1] - k2*y[1]; // pinene
  dydt[2] = k1 * y[1]; // dipentene
  dydt[3] = k2 * y[1] - 2*k4*y[3]*y[3] + 2*k5*y[5] -k3*y[3]; // alloocymene
  dydt[4] = k3*y[3]; // pyronene  
  dydt[5] = k4 * y[3]*y[3] - k5*y[5]; // dimer

  return dydt;                 }
}

data{
  int N_t;         // number of data
  int N_r;         // number of reactions
  int N_m;         // number of measured components
  array[N_t] vector[N_m] y;  // experimental data
  array[N_t] real ts;   // the time points at which to evaluate
  int N_p;        // the number of time points for prediction
  array[N_p] real t_pred; //time points for prediction
  real<lower=0> lkj_df; //parameter for lkj pior
}

parameters{
  real<lower=0> k1;
  real<lower=0> k2;
  real<lower=0> k3;
  real<lower=0> k4;
  real<lower=0> k5;
  vector[1] y0_initial;              // first y0 is estimated
  vector[N_m] sigma;    // each sigma is estimated
  cholesky_factor_corr[N_m] L_p; // cholesky factor for the correlation matrix
}
transformed parameters{
  vector[5] y0;
  y0[1]=y0_initial[1];
  y0[2]=0.0;
  y0[3]=0.0;
  y0[4]=0.0;
  y0[5]=0.0;
  array[N_t] vector[N_m] y_hat = ode_rk45(ode_pinene, 
                        y0, 0.0, ts, k1,k2,k3,k4,k5);

}
model{
  matrix[N_m,N_m] Sigma_chol = diag_pre_multiply(sigma,L_p);
  y0[1] ~ normal(100,10);
  k1 ~ lognormal(log(0.08),1);
  k2 ~ lognormal(log(0.04),1);
  k3 ~ normal(log(0.04),1);
  k4 ~ normal(log(0.06),1);
  k5 ~ normal(log(0.03),1);
  sigma ~ exponential(1);
  L_p ~ lkj_corr_cholesky(lkj_df);
  
        for (t in 1:N_t)
            y[t,]~multi_normal_cholesky(y_hat[t,], Sigma_chol[]);
}
 generated quantities{
 
  matrix[N_m, N_m] Omega = multiply_lower_tri_self_transpose(L_p); // correlation matrix
  matrix[N_m, N_m] Sigma = quad_form_diag(Omega, sigma);           // covariance matrix
  matrix[N_m, N_m] Sigma_chol = cholesky_decompose(Sigma);
  matrix[N_p,N_m] y_pred;
 
  // mean prediction for t_pred timepoints. This will give credible intervals for the mean prediction
  array[N_p] vector[N_m] y_hat2 = ode_rk45(ode_pinene, 
                        y0, 0.0, t_pred, k1,k2,k3,k4,k5);

  real cor_A_B = Omega[2,1];
  real cor_A_C = Omega[3,1];
  real cor_A_D = Omega[4,1];
  real cor_A_E = Omega[5,1];
  real cor_B_C = Omega[3,2];
  real cor_B_D = Omega[4,2];
  real cor_B_E = Omega[2,5];
  real cor_C_D = Omega[4,3];
  real cor_C_E = Omega[3,5];
  real cor_D_E = Omega[4,5];
  real<lower=0> sigma_prior;
  real k1_prior;
  real k2_prior;
  real k3_prior;
  real k4_prior;
  real k5_prior;
  real y0_prior;

  // prior simulations
  sigma_prior=exponential_rng(1);
  k1_prior = lognormal_rng(log(0.08),1);
  k2_prior = lognormal_rng(log(0.04),1);
  k3_prior = lognormal_rng(log(0.04),1);
  k4_prior = lognormal_rng(log(0.06),1);
  k5_prior = lognormal_rng(log(0.03),1);
  y0_prior = normal_rng(100,10);
  
  
   for (t in 1:N_p){
     y_pred[t,] = multi_normal_cholesky_rng(y_hat2[t,],Sigma_chol)';
   }

  }
