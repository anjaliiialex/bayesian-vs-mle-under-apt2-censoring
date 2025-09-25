
functions {
  real dagum_lpdf(real y, real a, real b, real p) {
    if (y <= 0) return negative_infinity();
    
    real log_y = log(y);
    real log_b = log(b);
    real y_over_b_neg_p = exp(-p * (log_y - log_b));
    
    return log(a) + log(p) - (p + 1) * log_y + p * log_b - (a + 1) * log1p(y_over_b_neg_p);
  }
  
  real dagum_lccdf(real y, real a, real b, real p) {
    if (y <= 0) return 0;
    
    real log_y = log(y);
    real log_b = log(b);
    real y_over_b_neg_p = exp(-p * (log_y - log_b));
    
    return -a * log1p(y_over_b_neg_p);
  }
}

data {
  int<lower=0> N_obs;
  vector<lower=0>[N_obs] y_obs;
  int<lower=0> R[N_obs];
  real<lower=0> p;
}

transformed data {
  real y_mean = mean(y_obs);
  real y_max = max(y_obs);
}

parameters {
  real<lower=0> a;
  real<lower=0> b;
}

model {
  // Informative but reasonable priors based on data
  // Center priors around plausible values
  a ~ normal(2, 1) T[0,];     // Truncated normal, positive
  b ~ normal(y_mean, y_mean) T[0,];  // Scale related to data scale
  
  // Likelihood
  for (i in 1:N_obs) {
    target += dagum_lpdf(y_obs[i] | a, b, p);
  }
  
  for (i in 1:N_obs) {
    if (R[i] > 0) {
      target += R[i] * dagum_lccdf(y_obs[i] | a, b, p);
    }
  }
}

generated quantities {
  real log_lik = 0;
  
  for (i in 1:N_obs) {
    log_lik += dagum_lpdf(y_obs[i] | a, b, p);
  }
  
  for (i in 1:N_obs) {
    if (R[i] > 0) {
      log_lik += R[i] * dagum_lccdf(y_obs[i] | a, b, p);
    }
  }
}

