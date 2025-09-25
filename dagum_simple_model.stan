// File: dagum_simple_model.stan
// A simplified diagnostic model for COMPLETE (non-censored) data -- CORRECTED

functions {
  real dagum_lpdf(real y, real a, real b, real p) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(a) + log(p) - (p + 1) * log(y) - (a + 1) * log(1 + (y/b)^(-p)) + p * log(b);
  }
}

data {
  int<lower=0> N;      
  vector<lower=0>[N] y; 
  real<lower=0> p;         
}

parameters {
  real log_a;
  real log_b;
}

transformed parameters {
  real<lower=0> a = exp(log_a);
  real<lower=0> b = exp(log_b);
}

model {
  // Priors on the unconstrained log-parameters
  log_a ~ normal(0.92, 0.5);
  log_b ~ normal(0.41, 0.5);

  // Stan doesn't have a built-in dagum distribution, so we must loop
  // through the data and add to the log-probability target manually.
  for (i in 1:N){
    target += dagum_lpdf(y[i] | a, b, p);
  }
}
