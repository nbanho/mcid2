# libraries
library(tidyverse)
library(rstan)
library(LaplacesDemon)
library(bayesplot)

# functions
softmax <- function(x) {
  exp(x) / sum(exp(x))
}

# seed
seed <- 12345
set.seed(seed)

# population and time
no_days <- 6 * 2 * 2
days <- 1:no_days
no_students <- 38
students <- 1:no_students

# input variables
class <- rep(c(0, 1), each = no_days/2)
aircleaner <- rep(c(0,1), times = 3, each = no_days/6)

# parameters
prob_test <- runif(no_days, 0.5, 1)
prob_virus <- c(0.5, 0.3, 0.2) # (SARS-CoV-2, Influenza, RSV)
no_virus <- length(prob_virus)
prob_pos <- 3 / no_students
beta0 <- logit(prob_pos * prob_virus)
beta1 <- c(1, 2, 1)
beta2 <- 0.0
beta3 <- 0.5

# initialize
susceptibles <- matrix(38, nrow = no_days + 1, ncol = no_virus)
y <- matrix(0, nrow = no_days, ncol = no_virus + 1)

# simulate data
mu <- numeric(no_virus + 1)
for (d in days) {
  tested_students <- sum(rbernoulli(no_students, prob_test[d]))
  mu[-1] <- beta0 - beta1 * aircleaner[d] + beta2 * class[d] + beta3 * (susceptibles[d] / 38)
  y[d,] <- rmultinom(1, tested_students, softmax(mu))
  susceptibles[d+1,] <- susceptibles[d,] - y[d,-1]
}

colSums(y)

# create stan data
sdl <- list(
  K = ncol(y),
  N = no_days,
  y = y,
  x1 = aircleaner,
  x2 = class,
  x3 = (susceptibles[-nrow(susceptibles),] / 38) / (2 * sd(susceptibles[-nrow(susceptibles),] / 38))
)


# model
mm <- stan(file = "tests/multi-logit.stan",
           data = sdl,
           seed = seed)

mm_sum <- summary(mm, pars = c("beta0", "beta1", "beta2", "beta3", "tau", "kappa"), probs = c(0.05, 0.95))
mm_sum$summary %>% data.frame() %>% mutate_if(is.numeric, round, 2)
mcmc_intervals(mm, pars = c("beta1[1]", "beta1[2]", "beta1[3]"))
mcmc_intervals(mm, pars = c("beta0[1]", "beta0[2]", "beta0[3]"))

