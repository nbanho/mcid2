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
  int C; // number of classes
  int samples[C*N,K]; // multinomial outcome for number of samples
  int aircleaner[C*N]; // binary indicator whether air cleaners were implemented
  int schoolclass[C*N]; // binary indicator for school class
  real susceptibles[C*N,K-1]; // number of susceptibles for each virus
}

transformed data {
  int K_c = K - 1; // #categories - 1
  int CN = C * N; // total number of observations across classes
}

parameters {
  vector[K_c] beta0; // intercept
  real beta1; // effect of air cleaners
  real beta2; // effect of class
  real beta3; // effect of susceptibles
}

transformed parameters {
  matrix[CN,K] mu;
  mu[,1] = rep_vector(0, CN);
  for (n in 1:CN) {
    for (k in 1:K_c) {
      mu[n,k+1] = beta0[k] + beta1 * aircleaner[n] + beta2 * schoolclass[n] + beta3 * susceptibles[n,k];
    }
  }
}

model {
  // priors
  beta0 ~ student_t(5., 0., 2.5);
  beta1 ~ student_t(5., 0., 5.);
  beta2 ~ student_t(5., 0., 5.);
  beta3 ~ student_t(5., 0., 5.);
  
  // likelihood
  for (n in 1:CN)
    samples[n] ~ multinomial_logit2_lpmf(mu[n]');
}