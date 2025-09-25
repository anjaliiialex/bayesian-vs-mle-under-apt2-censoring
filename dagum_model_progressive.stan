// File: dagum_model_progressive.stan
// FINAL VERSION with reparameterization to solve convergence issues

functions {
  real dagum_lpdf(real y, real a, real b, real p) {
    if (y <= 0) {
      return negative_infinity();
    }
    return log(a) + log(p) - (p + 1) * log(y) - (a + 1) * log(1 + (y/b)^(-p)) + p * log(b);
  }

  real dagum_lccdf(real y, real a, real b, real p) {
    if (y <= 0) {
      return 0; // log(1)
    }
    return -a * log(1 + (y/b)^(-p));
  }
}

data {
  int<lower=0> N_obs;      
  vector<lower=0>[N_obs] y_obs; 
  int<lower=0> R[N_obs];   
  real<lower=0> p;         
}

// Parameters are now defined on the unconstrained log scale.
parameters {
  real log_a; // log of parameter 'a'
  real log_b; // log of parameter 'b'
}

// We transform the log-parameters back to their natural scale
// so they can be used in the likelihood.
transformed parameters {
  real<lower=0> a = exp(log_a);
  real<lower=0> b = exp(log_b);
}

model {
  // --- PRIORS ---
  // Using the stronger priors to help with identifiability
  log_a ~ normal(0.92, 0.1);
  log_b ~ normal(0.41, 0.1);

  // --- LIKELIHOOD ---
  // (This part for the censored data remains the same)
  for (i in 1:N_obs) {
    target += dagum_lpdf(y_obs[i] | a, b, p);
  }

  for (i in 1:N_obs) {
    if (R[i] > 0) {
      target += R[i] * dagum_lccdf(y_obs[i] | a, b, p);
    }
  }
}

