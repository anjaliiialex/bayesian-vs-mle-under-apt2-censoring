// File: dagum_model_progressive_fixed.stan
// Fixed version with proper uninformative priors and numerical stability

functions {
  // More numerically stable Dagum PDF
  real dagum_lpdf(real y, real a, real b, real p) {
    if (y <= 0) {
      return negative_infinity();
    }
    
    real log_y = log(y);
    real log_b = log(b);
    real log_ratio = -p * (log_y - log_b);  // log((y/b)^(-p))
    
    // Avoid overflow in log1p_exp for large values
    real log_1_plus_ratio;
    if (log_ratio > 20) {
      log_1_plus_ratio = log_ratio;  // log(1 + exp(x)) ≈ x for large x
    } else {
      log_1_plus_ratio = log1p_exp(log_ratio);  // log(1 + exp(x))
    }
    
    return log(a) + log(p) - (p + 1) * log_y - (a + 1) * log_1_plus_ratio + p * log_b;
  }
  
  // More numerically stable survival function
  real dagum_lccdf(real y, real a, real b, real p) {
    if (y <= 0) {
      return 0; // log(1)
    }
    
    real log_y = log(y);
    real log_b = log(b);
    real log_ratio = -p * (log_y - log_b);  // log((y/b)^(-p))
    
    // Numerical stability for log(1 + x)
    real log_1_plus_ratio;
    if (log_ratio > 20) {
      log_1_plus_ratio = log_ratio;
    } else {
      log_1_plus_ratio = log1p_exp(log_ratio);
    }
    
    return -a * log_1_plus_ratio;
  }
}

data {
  int<lower=0> N_obs;      
  vector<lower=0>[N_obs] y_obs; 
  int<lower=0> R[N_obs];   
  real<lower=0> p;         
}

parameters {
  real log_a; // log of parameter 'a'
  real log_b; // log of parameter 'b'
}

transformed parameters {
  real<lower=0> a = exp(log_a);
  real<lower=0> b = exp(log_b);
}

model {
  // --- PROPER UNINFORMATIVE PRIORS ---
  // These are much more reasonable and don't give unfair advantage
  log_a ~ normal(0, 1.5);  // Allows a wide range of values for 'a'
  log_b ~ normal(0, 1.5);  // Allows a wide range of values for 'b'
  
  // Alternative: Even more uninformative
  // log_a ~ normal(0, 2.5);
  // log_b ~ normal(0, 2.5);
  
  // --- LIKELIHOOD ---
  // Likelihood for observed failures
  for (i in 1:N_obs) {
    target += dagum_lpdf(y_obs[i] | a, b, p);
  }
  
  // Likelihood for progressively censored units
  for (i in 1:N_obs) {
    if (R[i] > 0) {
      target += R[i] * dagum_lccdf(y_obs[i] | a, b, p);
    }
  }
}
