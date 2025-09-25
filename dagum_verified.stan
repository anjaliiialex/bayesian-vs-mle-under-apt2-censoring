
functions {
  real dagum_lpdf(real x, real a, real b, real p) {
    if (x <= 0) return negative_infinity();
    
    real log_pdf = log(a) + log(p) - (p + 1) * log(x) + p * log(b) - (a + 1) * log1p(pow(x/b, -p));
    return log_pdf;
  }
  
  real dagum_lccdf(real x, real a, real b, real p) {
    if (x <= 0) return 0;
    
    return -a * log1p(pow(x/b, -p));
  }
}

data {
  int<lower=0> N_obs;
  vector<lower=0>[N_obs] y_obs;
  int<lower=0> R[N_obs];
  real<lower=0> p;
}

parameters {
  real<lower=0> a;
  real<lower=0> b;
}

model {
  // Weakly informative priors
  a ~ gamma(2, 1);  // mean = 2, variance = 2
  b ~ gamma(2, 1);  // mean = 2, variance = 2
  
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

