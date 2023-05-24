functions {
  // generated with brms 2.18.0
  /* multinomial-logit log-PMF
   * Args:
   *   y: array of integer response values
   *   mu: vector of category logit probabilities
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real multinomial_logit2_lpmf(int[] y, vector mu) {
     return multinomial_lpmf(y | softmax(mu));
   }
}

data {
  int K; // number of categories
  int N; // number of observations
  int y[N,K]; // multinomial outcome
  int x1[N]; // aircleaner
  int x2[N]; // class
  real x3[N,K-1]; // susceptibles
}

transformed data {
  int K_c = K - 1; // #categories - 1
}

parameters {
  vector[K_c] beta0; // intercept
  vector[K_c] beta1; // effect of air cleaners
  real kappa; // average effect of air cleaners across pathogens
  real<lower=0> tau; // variation in effect of air cleaners
  real beta2; // effect of class
  real beta3; // effect of susceptibles
}

transformed parameters {
  matrix[N,K] mu;
  mu[,1] = rep_vector(0, N);
  for (n in 1:N) {
    for (k in 1:K_c) {
      mu[n,k+1] = beta0[k] + beta1[k] * x1[n] + beta2 * x2[n] + beta3 * x3[n,k];
    }
  }
}

model {
  // priors
  beta0 ~ student_t(5, 0, 2.5);
  kappa ~ student_t(5, 0, 5);
  tau ~ student_t(5,0,1);
  beta1 ~ normal(kappa, tau);
  beta2 ~ student_t(5, 0, 5);
  beta3 ~ student_t(5, 0, 5);
  
  // likelihood
  for (n in 1:N)
    y[n] ~ multinomial_logit2_lpmf(mu[n]');
}