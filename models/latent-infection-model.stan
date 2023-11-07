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
  int<lower=1> SS; // number of days from last Saturday to study start
  int<lower=1> W; // number of weeks
  int<lower=1> C; // number of classes
  int<lower=1> K; // number of viruses
  int cases[D,C]; // number of new cases per class
  vector<lower=0>[7] totCasesWeekday; // number of cases by weekday
  int<lower=0,upper=1> weekend[D+S]; // weekend indicator
  int<lower=0,upper=1> vacation[D+S]; // vacation indicator
  int<lower=0,upper=1> airCleaner[D+S,C]; // air filter indicator per class
  matrix[D,C] stud; // proportion of students in class per class
  matrix[D,C] vent; // air change rate (ventilation) per class
  int<lower=0> Mvent; // number of missing values in aer
  int<lower=0,upper=Mvent> missingVent[D,C]; // missing aer change rates indicator per class
  vector[D+S] cov; // proportion of positive covid-19 tests in Solothurn
  vector[D+S] ili; // number of consultations regarding influenza-like illnesses in Solothurn
  real p_in_mu_m[K]; // prior: location hyperparameter m in mu^p_IN ~ Normal(m, s)
  real pos_samples[D,K,C]; // weekly number of positive samples by virus
  real p_in_mu_s[K]; // prior: scale hyperparameter s in mu^p_IN ~ Normal(m, s)
  real p_in_sigma_m[K]; // prior: location hyperparameter m in sigma^p_IN ~ Normal(m, s)
  real p_in_sigma_s[K]; /// prior: scale hyperparameter m in sigma^p_IN ~ Normal(m, s)
}

transformed data {
  int week_index[W,2];
  real p_weights[D,K,C];
  for (w in 1:W) {
    week_index[w,1] = (w - 1) * 7 + 1;
    week_index[w,2] = w * 7;
  }
  for (c in 1:C) {
    for (d in 1:D) {
      for (k in 1:K) {
        p_weights[d,k,c] = pos_samples[d,k,c] / sum(pos_samples[d,,c]);
      }
    }
  }
}


parameters {
  simplex[7] vega; // weights of weekday cases
  real<lower=0> invphi;
  real alpha; // model intercept
  real omega; // weekend effect
  vector[6] beta; // aircleaner effect
  vector[Mvent] ventNa; // estimates for missing aer
  vector[K] mu_p_in; // log mean in p_IN ~ Lognormal(mu, sigma)
  vector<lower=0>[K] sigma_p_in; // log standard deviation in p_IN ~ Lognormal(mu, sigma)
  vector<lower=0>[C] I0; // number of infections at the first seeding day
}


transformed parameters {
  real<lower=0> phi = inv_square(invphi); // over-dispersion parameter
  matrix<lower=0>[D+SS,C] muCases; // expected number of new cases
  matrix[D+S,C] logInfections; // log of the expected number of new infections
  matrix<lower=0>[D+S,C] infections; // log of the expected number of new infections
  matrix<lower=0>[D+S,C] cumInfections; // cumulative number of infections
  matrix<lower=0>[D+S,C] contagious; // number of contagious/infectious students
  matrix<lower=0>[D+S,K] p_in; // probability distribution for incubation periods
  
  // compute discretized p_IN distributions
  for (k in 1:K) {
    for (d in 1:(D+S)) { 
      p_in[d,k] = diff_lnorm(D+S-d, mu_p_in[k], sigma_p_in[k]);
    }
  }
  
  for (c in 1:C) {
    
    // seeding phase
    contagious[1,c] = 0;
    cumInfections[1,c] = 0;
    logInfections[1,c] = log(I0[c]);
    infections[1,c] = I0[c];
    for (d in 2:S) {
      // contagious
      contagious[d,c] = sum(infections[max(1,(d-7)):(d-1),c]);
      // compute cumulative infections
      cumInfections[d,c] = sum(infections[1:(d-1),c]);
      // estimate log number of infections
      logInfections[d,c] = log(contagious[d,c] / cumInfections[d,c]) + alpha + omega * fmax(weekend[d], vacation[d]) + beta[2] * (c-1) + beta[5] * cov[d] +  beta[6] * ili[d];
      // ccompute number of infections
      infections[d,c] = exp(logInfections[d,c]);
    }
    
    // expected number of cases SS days before study start
    for (d in 1:SS) {
      muCases[d,c] = dot_product(infections[1:(S-SS+d),c], block(p_in, D+SS-d+1, 1, S+d-2, K) * rep_vector(1. / (1. * K), K));
    }
  
    // compute infections and expected number of cases
    for (d in (S+1):(D+S)) {
      // contagious 
      contagious[d,c] = sum(infections[(d-7):(d-1),c]);
      // compute cumulative infections
      cumInfections[d,c] = sum(infections[1:(d-1),c]);
      // estimate log number of infections on every day
      logInfections[d,c] = log(contagious[d,c] / cumInfections[d,c]) + alpha + omega * fmax(weekend[d], vacation[d]) + beta[2] * (c-1) + beta[5] * cov[d] +  beta[6] * ili[d];
      // log number of infections during school days
      if (fmax(weekend[d], vacation[d]) == 0) {
        // considering missing data points in vent
        if (missingVent[d-S,c] > 0) { 
          logInfections[d,c] = logInfections[d,c] + beta[1] * airCleaner[d,c] + beta[3] * stud[d-S,c] + beta[4] * ventNa[missingVent[d-S,c]];
        } else {
          logInfections[d,c] = logInfections[d,c] + beta[1] * airCleaner[d,c] + beta[3] * stud[d-S,c] + beta[4] * vent[d-S,c];
        }
      }
      // compute expected number of new cases
      infections[d,c] = exp(logInfections[d,c]);
      muCases[d-S+SS,c] = dot_product(infections[1:d,c], block(p_in, S+D-d+1, 1, d, K) * to_vector(p_weights[d-S,,c]));
    }
  
    // re-weight expected cases to consider weekday effects
    for (w in 1:W) {
      muCases[week_index[w,1]:week_index[w,2],c] = re_weight(muCases[week_index[w,1]:week_index[w,2],c], vega);
    }
    
  }
}


model {
  // priors
  vega ~ dirichlet(totCasesWeekday);
  invphi ~ normal(0., 1.);
  alpha ~ student_t(5., 0., 10.);
  omega ~ normal(log(1.1), .05);
  mu_p_in ~ normal(p_in_mu_m, p_in_mu_s);
  sigma_p_in ~ normal(p_in_sigma_m, p_in_sigma_s);
  beta ~ student_t(5., 0., 2.5);
  ventNa ~ normal(0., 1.);
  I0 ~ exponential(1.);
  
  // likelihood 
  for (c in 1:C) {
    for (d in 1:D) {
      if (vacation[S+d] == 0) {
        cases[d,c] ~ neg_binomial_2(muCases[d+SS,c], phi);
        //cases[d,c] ~ poisson(muCases[d+SS,c]);
      }
    }
  }
}

generated quantities {
  real contagious_ac[D+S,C,3]; // number of contagious individuals
  real cum_infections_ac[D+S,C,3]; // cum infections depending by installation of air cleaners
  real log_infections_ac[D+S,C,3]; // log infections depending by installation of air cleaners
  real infections_ac[D+S,C,3]; // infections depending by installation of air cleaners
  real mu_cases_ac[D+SS,C,3]; // expected cases depending by installation of air cleaners

  for (a in 1:3) {
    for (c in 1:C) {

      // seeding phase
      contagious_ac[1,c,a] = 0;
      cum_infections_ac[1,c,a] = 0;
      log_infections_ac[1,c,a] = log(I0[c]);
      infections_ac[1,c,a] = I0[c];
      for (d in 2:S) {
        // contagious
        contagious_ac[d,c,a] = sum(infections_ac[max(1,(d-7)):(d-1),c,a]);
        // compute cumulative infections
        cum_infections_ac[d,c,a] = sum(infections_ac[1:(d-1),c,a]);
        // estimate log number of infections
        log_infections_ac[d,c,a] = log(contagious_ac[d,c,a] / cum_infections_ac[d,c,a]) + alpha + omega * fmax(weekend[d], vacation[d]) + beta[2] * (c-1) + beta[5] * cov[d] +  beta[6] * ili[d];
        // ccompute number of infections
        infections_ac[d,c,a] = exp(log_infections_ac[d,c,a]);
      }

      // expected number of cases SS days before study start
      for (d in 1:SS) {
        mu_cases_ac[d,c,a] = dot_product(to_vector(infections_ac[1:(S-SS+d),c,a]), block(p_in, D+SS-d+1, 1, S+d-2, K) * rep_vector(1. / (1. * K), K));//tail(to_vector(p_in_eff[d,,c]), S-SS+d));
      }

      // compute infections and expected number of cases
      for (d in (S+1):(D+S)) {
        contagious_ac[d,c,a] = sum(infections_ac[(d-7):(d-1),c,a]);
        // compute cumulative infections
        cum_infections_ac[d,c,a] = sum(infections_ac[1:(d-1),c,a]);
        // estimate log number of infections on every day
        log_infections_ac[d,c,a] = log(contagious_ac[d,c,a] / cum_infections_ac[d,c,a]) + alpha + omega * fmax(weekend[d], vacation[d]) + beta[2] * (c-1) + beta[5] * cov[d] +  beta[6] * ili[d];
        // log number of infections during school days
        if (fmax(weekend[d], vacation[d]) == 0) {
          if (a == 1) { // air cleaners installed according to study design
            log_infections_ac[d,c,a] = log_infections_ac[d,c,a] + beta[1] * airCleaner[d,c];
          }
          if (a == 3) { // air cleaners fully installed over the study
            log_infections_ac[d,c,a] = log_infections_ac[d,c,a] + beta[1];
          }
          if (missingVent[d-S,c] > 0) { // considering missing data points in vent
            log_infections_ac[d,c,a] = log_infections_ac[d,c,a] + beta[3] * stud[d-S,c] + beta[4] * ventNa[missingVent[d-S,c]];
          } else {
            log_infections_ac[d,c,a] = log_infections_ac[d,c,a] + beta[3] * stud[d-S,c] + beta[4] * vent[d-S,c];
          }
        }
        // compute expected number of new cases
        infections_ac[d,c,a] = exp(log_infections_ac[d,c,a]);
        mu_cases_ac[d-S+SS,c,a] = dot_product(to_vector(infections_ac[1:d,c,a]), block(p_in, S+D-d+1, 1, d, K) * to_vector(p_weights[d-S,,c]));//tail(to_vector(p_in_eff[d,,c]), d));
      }

      // re-weight expected cases to consider weekday effects
      for (w in 1:W) {
        mu_cases_ac[week_index[w,1]:week_index[w,2],c,a] = to_array_1d(re_weight(to_vector(mu_cases_ac[week_index[w,1]:week_index[w,2],c,a]), vega));
      }
    }
  }
}
