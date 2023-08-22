data {
  int<lower=1> Da; // number of days in class A
  int<lower=1> Db; // number of days in class B
  int coughsA[Da]; // number of coughs in class A
  int coughsB[Db]; // number of coughs in class B
  int tA[Da]; // duration of cough detection in class A
  int tB[Db]; // duration of cough detection in in class B
  int<lower=0,upper=1> airCleanerA[Da]; // air filter indicator in class A
  int<lower=0,upper=1> airCleanerB[Db]; // air filter indicator in class B
  real nA[Da]; // number of students persent in class A (standardized)
  real nB[Db]; // number of students persent in class A (standardized)
  real aerA[Da]; // air change rate in class A (standardized)
  real aerB[Db]; // air change rate in class B (standardized)
  real casesA[Da]; // cumulative cases in class A (standardized)
  real casesB[Db]; // cumulative cases in class A (standardized)
  real weekdayA[Da,4]; // weekday in class A (standardized)
  real weekdayB[Db,4]; // weekday in class A (standardized)
  real s_x[9]; // prior scale adjustments for coefficients (aircleaner, class, Tue, Wed, Thu, Fri, n, aer, cases)
}


parameters {
  real<lower=0> invphi; // reciprocal over-dispersion parameter
  real beta0; // model intercept
  real beta1; // aircleaner effect
  real beta2; // class effect
  vector[4] beta3; // weekday effects 
  real beta4; // students effect
  real beta5; // ventilation / air quality effect
  real beta6; // cumulative cases effect
}


transformed parameters {
  real<lower=0> phi = inv_square(invphi); // over-dispersion parameter
  vector[Da] log_muA; // expected coughs in class A
  vector[Db] log_muB; // expected coughs in class B
  
  for (i in 1:Da) {
    log_muA[i] = log(tA[i]) + beta0 + beta1 * airCleanerA[i] + beta3[1] * weekdayA[i,1] + beta3[2] * weekdayA[i,2] + beta3[3] * weekdayA[i,3] + beta3[4] * weekdayA[i,4] + beta4 * nA[i] + beta5 * aerA[i] + beta6 * casesA[i];
  }
  for (j in 1:Db) {
    log_muB[j] = log(tB[j]) + beta0 + beta1 * airCleanerB[j] + beta2 + beta3[1] * weekdayB[j,1] + beta3[2] * weekdayB[j,2] + beta3[3] * weekdayB[j,3] + beta3[4] * weekdayB[j,4] + beta4 * nB[j] + beta5 * aerB[j] + beta6 * casesB[j];
  }
}


model {
  // priors
  invphi ~ normal(0., 1.);
  beta0 ~ student_t(5., 0, 2.5);
  beta1 ~ student_t(5., 0, 2.5 / s_x[1]);
  beta2 ~ student_t(5., 0, 2.5 / s_x[2]);
  for (wd in 1:4) {
    beta3[wd] ~ student_t(5., 0, 2.5 / s_x[2+wd]);
  }
  beta4 ~ student_t(5., 0, 2.5 / s_x[7]);
  beta5 ~ student_t(5., 0, 2.5 / s_x[8]);
  beta6 ~ student_t(5., 0, 2.5 / s_x[9]);
  
  // likelihood 
  for (i in 1:Da) {
    coughsA[i] ~ neg_binomial_2_log(log_muA[i], phi);
  }
  for (j in 1:Db) {
    coughsB[j] ~ neg_binomial_2_log(log_muB[j], phi);
  }
}
