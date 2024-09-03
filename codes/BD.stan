// Stan code for the Box-Draper model with estimation of the decomposed covar matrix using cholesky factorization
// This code adds prior simularion in the generated quantities

data{
  int<lower=1> N_t;                  // number of data
  int<lower=1> N_r;                  // number of responses
  int<lower=1> N_p;                 // number of time points for prediction  
  array[N_t] vector[N_r] y;   //three measured components
  vector[N_t] ts;      // the time points at which to evaluate
  array[N_p] real t_pred; //time points for prediction
  real<lower=0> lkj_df;    //parameter for lkj prior
}

parameters{
  real<lower=0> k1;
  real<lower=0> k2;
  vector[1] y0_initial;       // first y0 is estimated
  vector[N_r] sigma;          // standard deviation for each response
  cholesky_factor_corr[N_r] L_p; // cholesky factor for the correlation matrix
}

transformed parameters{
    vector[3] y0;
    y0[1] = y0_initial[1];
    y0[2] = 0.0;
    y0[3] = 0.0;
 
    array[N_t] vector[N_r] y_hat;
    
    for (t in 1:N_t){
       y_hat[t,1] = y0[1]*exp(-k1*ts[t]);
       y_hat[t,2] = y0[1]*k1/(k2-k1)*(exp(-k1*ts[t])-exp(-k2*ts[t]));
       y_hat[t,3] = y0[1]*(1+(k1*exp(-k2*ts[t])-k2*exp(-k1*ts[t]))/(k2-k1));
   }

}

model{
  
  matrix[N_r,N_r] Sigma_chol = diag_pre_multiply(sigma,L_p);
  
  y0[1] ~ normal(1,0.1);
  k1 ~ lognormal(log(0.2), 1);
  k2 ~ lognormal(log(0.5),1);
  sigma ~ exponential(1);
  L_p ~ lkj_corr_cholesky(lkj_df);

  y ~ multi_normal_cholesky(y_hat, Sigma_chol);
}

generated quantities{
  matrix[N_r, N_r] Omega = multiply_lower_tri_self_transpose(L_p); // correlation matrix
  matrix[N_r, N_r] Sigma = quad_form_diag(Omega, sigma);           // covariance matrix
  
  real cor_A_B = Omega[2,1];
  real cor_A_C = Omega[3,1];
  real cor_B_C = Omega[3,2];
  matrix[N_p, N_r] y_hat2;  // mean prediction for t_pred timepoints. This will give credible intervals for the mean prediction
  matrix[N_p,N_r] y_pred;
  matrix[N_t, N_r] y_rep;   // Generate posterior predictive checks 
  real<lower=0> sigma_prior;
  real k1_prior;
  real k2_prior;
  real A0_prior;

  // prior simulations
  sigma_prior=exponential_rng(1);
  k1_prior = lognormal_rng(log(0.2),1);
  k2_prior = lognormal_rng(log(0.5),1);
  A0_prior = normal_rng(1,0.1);

  {
    
    matrix[N_r, N_r] Sigma_chol = cholesky_decompose(Sigma);

   for (t in 1:N_p){
       y_hat2[t,1] = y0[1]*exp(-k1*t_pred[t]);
       y_hat2[t,2] = y0[1]*k1/(k2-k1)*(exp(-k1*t_pred[t])-exp(-k2*t_pred[t]));
       y_hat2[t,3] = y0[1]*(1+(k1*exp(-k2*t_pred[t])-k2*exp(-k1*t_pred[t]))/(k2-k1));
   }
   
   for (t in 1:N_p){
     y_pred[t,] = multi_normal_cholesky_rng(y_hat2[t,],Sigma_chol)';
   }

      for(t in 1:N_t){
      y_rep[t,] = multi_normal_cholesky_rng(y_hat[t,], Sigma_chol)';

    }
  }
}
