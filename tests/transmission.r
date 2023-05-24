# libraries
library(tidyverse)
library(rstanarm)
library(LaplacesDemon)

# functions
source("utils/epi.r")

re_weight <- function(x, p) {
  x * p * (sum(x) / (sum(x * p)))
}

# seed
seed <- 12345
set.seed(seed)

# population and time
weeks <- 7
S <- 40
D <- weeks * 7
weekday <- c(c("Sat", "Sun"),
             rep(c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"), times = 6),
             c("Mon", "Tue", "Wed", "Thu", "Fri"))
days <- 1:D
week <- c(c(0,0),
          rep(1:6, each = 7),
          rep(7, 5))

# input variables
aircleaner <- ifelse(!(weekday %in% c("Sat", "Sun")) & week %in% c(3, 5, 6), 1, 0)
weekend <- ifelse(weekday %in% c("Sat", "Sun"), 1, 0)
vacation <- ifelse(week == 4, 1, 0)
noschool <- pmax(vacation, weekend)

# parameters
t0 <- 7
I0 <- rexp(t0, rate = t0)
alpha <- -2
omega <- -0.15
beta1 <- -1
pIN <- vp(D + t0, 0, p_in)
w <- c(0.025, 0.025, 0.35, 0.15, 0.15, 0.15, 0.15)


# simulate data 
muC <- numeric(D)
I <- numeric(D)
for (d in 1:D) {
  if (d == 1) {
    N <- sum(I0)
  } else {
    N <- sum(c(I0, I[1:(d-1)]))
  }
  I[d] <- N * exp(alpha + omega * noschool[d] + beta1 * aircleaner[d])
  muC[d] <- sum(rev(pIN[1:(d+t0)]) * c(I0, I[1:d]))
}

# cases without shift
C <- sapply(muC, rpois, n = 1)

# cases with shift
muC_weekly <- split(muC, ceiling(seq_along(muC)/7))
muC_tilde <- unlist(lapply(muC_weekly, re_weight, p = w))
C_tilde <- sapply(muC_tilde, rpois, n = 1)

# inspect
df <- data.frame(
  C = C,
  C_tilde = C_tilde,
  aircleaner = aircleaner,
  day = 1:D,
  weekday = WD
)

df %>%
  group_by(weekday) %>%
  summarize(nC = sum(C),
            nCtilde = sum(C_tilde))

nCtildeWeekly <- df %>% 
  group_by(weekday) %>% 
  summarize(n = sum(C_tilde)) %>% 
  ungroup() %>% 
  mutate(weekday = factor(weekday, levels = c("Sat", "Sun", "Mon", "Tue", "Wed", "Thu", "Fri"))) %>%
  arrange(weekday) %>%
  dplyr::select(n) %>% 
  unlist

# stan data
sdl <- list(
  S = t0,
  D = D,
  W = weeks,
  Cases = C_tilde,
  totCasesWeekly = nCtildeWeekly + 1,
  NoSchool = noschool,
  AirCleaner = aircleaner,
  p_in_mu_m = p_in_mu_m,
  p_in_mu_s = p_in_mu_s,
  p_in_sigma_m = p_in_sigma_m,
  p_in_sigma_s = p_in_sigma_s
)

# model
lnb <- stan(file = "tests/latent-negbinom.stan", 
            data = sdl,
            seed = seed)

lnb_sum <- summary(lnb, pars = c("alpha", "omega", "beta1", "I0", "vega", "mu_p_in", "sigma_p_in"), probs = c(0.05, 0.95))
lnb_sum$summary %>% data.frame() %>% mutate_if(is.numeric, round, 2)

mu_cases <- summary(lnb, pars = "muCases")$summary[,"mean"]

result <- data.frame(
  true_cases = C_tilde,
  expected_cases = mu_cases,
  weekday = weekday
)

result %>% 
  group_by(weekday) %>% 
  summarize(n_true = sum(true_cases),
            n_esti = sum(expected_cases)) %>% 
  ungroup() %>% 
  mutate(weekday = factor(weekday, levels = c("Sat", "Sun", "Mon", "Tue", "Wed", "Thu", "Fri"))) %>%
  arrange(weekday) 
  
