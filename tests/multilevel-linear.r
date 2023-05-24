library(tidyverse)
library(rstanarm)
library(LaplacesDemon)

aircleaner <- rep(c(0,1), each = 6)
pathogen <- rep(c("SARS-CoV-2", "Influenza", "RSV"), each = 5)
dat <- expand_grid(aircleaner = aircleaner, pathogen = factor(pathogen)) 
binary_pathogen <- model.matrix(~ pathogen - 1, data = dat)
dat <- cbind(dat, binary_pathogen) %>%
  mutate(mu_y = 2 - 0.5 * aircleaner - .25 * pathogenInfluenza + 0.5 * `pathogenSARS-CoV-2`,
         y = sapply(mu_y, function(m) rtrunc(n = 1, spec = "norm", a = 0, mean = m, sd = .5))) 


mlm <- stan_glmer(y ~ 1 + aircleaner + (1 + aircleaner | pathogen), data = dat)

summary(mlm)
