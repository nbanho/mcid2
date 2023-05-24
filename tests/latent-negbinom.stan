functions {
  // Discretization of the cumulative lognormal distribution function
  real diff_lnorm(real x, real mu, real sigma) {
    if (x == 0) {
      return lognormal_cdf(0.5, mu, sigma);
    } else {
      return lognormal_cdf(x+0.5, mu, sigma) - lognormal_cdf(x-0.5, mu, sigma);
    }
  }
  
  // re-weight cases
  vector re_weight(vector x, vector p) {
    int K = num_elements(x);
    vector[K] y;
    for (k in 1:7) {
      y[k] = x[k] * p[k] * (sum(x) / dot_product(x, p));
    }
    return(y);
  }
}


data {
  int<lower=1> D; // number of days
  int<lower=1> S; // number of seeding days
  int<lower=1> W; // number of weeks
  int Cases[D]; // number of cases
  vector[7] totCasesWeekly; // number of cases by weekday
  int NoSchool[D]; // weekend and vacation indicator
  int AirCleaner[D]; // air filter indicator
  real p_in_mu_m; // prior: location hyperparameter m in mu^p_IN ~ Normal(m, s)
  real p_in_mu_s; // prior: scale hyperparameter s in mu^p_IN ~ Normal(m, s)
  real p_in_sigma_m; // prior: location hyperparameter m in sigma^p_IN ~ Normal(m, s)
  real p_in_sigma_s; /// prior: scale hyperparameter m in sigma^p_IN ~ Normal(m, s)
}

transformed data {
  int week_index[W,2];
  for (w in 1:W) {
    week_index[w,1] = (w - 1) * 7 + 1;
    week_index[w,2] = w * 7;
  }
}


parameters {
  real alpha; // model intercept
  real omega; // weekend effect
  real beta1; // aircleaner effect
  simplex[7] vega; // weights of weekday cases
  vector<lower=0>[S] I0; // cases in the seeding period
  real mu_p_in; // log mean in p_IN ~ Lognormal(mu, sigma)
  real<lower=0> sigma_p_in; // log standard deviation in p_IN ~ Lognormal(mu, sigma)
}


transformed parameters {
  vector<lower=0>[D] muCases; // expected number of new cases
  vector<lower=0>[D] Infections; // expected number of new infections
  vector<lower=0>[D] cumInfections; // cumulative number of infections
  vector<lower=0>[D+S] p_in; // probability distribution for incubation period
  
  // Compute discretized p_IN distribution
  for (k in 1:(D+S)) { 
    p_in[k] = diff_lnorm(D+S-k, mu_p_in, sigma_p_in);
  }
  
  // compute infections and expected number of cases
  for (d in 1:D) {
    cumInfections[d] = sum(append_row(I0, Infections[1:(d-1)]));
    Infections[d] = cumInfections[d] * exp(alpha + omega * NoSchool[d] + beta1 * AirCleaner[d]);
    muCases[d] = dot_product(append_row(I0, Infections[1:d]), tail(p_in, S+d));
  }
  
  // re-weight expected cases
  for (w in 1:W) {
    muCases[week_index[w,1]:week_index[w,2]] = re_weight(muCases[week_index[w,1]:week_index[w,2]], vega);
  }
}


model {
  // priors
  alpha ~ student_t(5., 0., 10.);
  omega ~ normal(-0.15, .15);
  mu_p_in ~ normal(p_in_mu_m, p_in_mu_s);
  sigma_p_in ~ normal(p_in_sigma_m, p_in_sigma_s);
  beta1 ~ student_t(5., 0., 2.5);
  I0 ~ exponential(S / 1.);
  vega ~ dirichlet(totCasesWeekly);

  // likelihood 
  Cases ~ poisson(muCases);
}