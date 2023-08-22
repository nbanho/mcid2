data {
  int<lower=0> N; // number of observations
  int<lower=1> V; // number of viruses
  int<lower=0> coughs[N]; // number of coughs
  int<lower=0> TT[N]; // number of minutes
  int<lower=0> viruses[N,V]; // number of positive tests for each virus
  real s_v; // prior scale for theta
}

parameters {
  real<lower=0> invphi; // inverse of over-dispersion parameter
  real theta; // average effect of a positive test
  real<lower=0> tau; // variation in effect of a positive test
  vector[V] theta_v; // virus-specific effect of a positive test
  real beta0; // model intercept
}

transformed parameters {
  real<lower=0> phi = inv_square(invphi); // over-dispersion parameter
  vector[N] log_mu; // mean parameters of negative binomial model
  
  for (i in 1:N) {
    log_mu[i] = log(TT[i]) + beta0 + dot_product(theta_v, to_vector(viruses[i,1:V]));
  }
}

model {
  // priors
  invphi ~ normal(0., 1.);
  beta0 ~ student_t(5., 0, 2.5);
  
  tau ~ normal(0., 1.);
  theta ~ student_t(5., 0., 2.5 / s_v);
  theta_v ~ normal(theta, tau);
  
  // likelihood
  for (i in 1:N) {
    coughs[i] ~ neg_binomial_2_log(log_mu[i], phi);
  }
}

