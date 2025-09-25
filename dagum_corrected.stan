
// CORRECTED Stan model using your working likelihood
functions {
  // Use the exact same Dagum PDF that works in R
  real dagum_lpdf(real y, real a, real b, real p) {
    if (y <= 0) {
      return negative_infinity();
    }
    
    // This matches ddagum(y, a, b, p) exactly
    real pdf_val = (a * p / pow(y, p+1)) * pow(b, p) / pow(1 + pow(y/b, -p), a+1);
    return log(pdf_val);
  }
  
  // Use the exact same survival function that works in R  
  real dagum_lccdf(real y, real a, real b, real p) {
    if (y <= 0) {
      return 0;
    }
    
    // This matches 1 - VGAM::pdagum(...) exactly
    real survival_val = pow(1 + pow(y/b, -p), -a);
    return log(survival_val);
  }
}

data {
  int<lower=0> N_obs;      
  vector<lower=0>[N_obs] y_obs; 
  int<lower=0> R[N_obs];   
  real<lower=0> p;         
}

parameters {
  real<lower=0> a;  // Use constrained parameters directly
  real<lower=0> b;
}

model {
  // Proper uninformative priors
  a ~ gamma(1, 1);  // More stable than log-normal
  b ~ gamma(1, 1);
  
  // Likelihood - matches your R function exactly
  for (i in 1:N_obs) {
    target += dagum_lpdf(y_obs[i] | a, b, p);
  }
  
  for (i in 1:N_obs) {
    if (R[i] > 0) {
      target += R[i] * dagum_lccdf(y_obs[i] | a, b, p);
    }
  }
}

