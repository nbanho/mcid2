#### Library ####

library(tidyverse)
library(reshape2)
library(lubridate)


#### Data ####

files_A <- list.files("data-raw/palas/E3a_komplett_aqGuard_15988", full.names = T)
files_B <- list.files("data-raw/palas/E3b_komplett_aqGuard_15973", full.names = T)
files_A <- files_A[grepl("txt", files_A)]
files_B <- files_B[grepl("txt", files_B)]

read_palas <- function(file) {
  df <- read.table(file, sep ="\t")
  colnames(df) <- df[1,]
  df <- df[-1,c(1, 3:10, 19, 21)]
  colnames(df) <- c("datetime", "PM1", "PM25", "PM4", "PM10", "PMtot", "CO2", "VOC", "Cn", "Temp", "rH")
  
  df <- df %>%
    mutate(across(-datetime, as.numeric),
           date = as.Date(datetime, format = "%d.%m.%Y"),
           datetime = as.POSIXct(datetime, format = "%d.%m.%Y %H:%M:%S"),
           time = format(datetime, format = "%H:%M:%S")) %>%
    select(date, time, datetime, everything())
  
  return(df)
}


df <- rbind(data.frame(class = "A", file = files_A), 
            data.frame(class = "B", file = files_B)) %>%
  mutate(data = lapply(file, read_palas)) %>%
  unnest(data)


#### Preprocessing ####

# remove weekends and vacation
df <- df %>%
  filter(!(weekdays(date) %in% c("Saturday", "Sunday")),
         !(date %in% seq.Date(as.Date("2023-02-06"), as.Date("2023-02-11"), by = "1 day")),
         date >= as.Date("2023-01-16"),
         date <= as.Date("2023-03-11")) 

# move time in class B (E3b) on 2023-01-16 (first weekend, Monday) by one hour backwards (see data-collection-issues.docx)
df$datetime[df$class=="B"&df$date==as.Date("2023-01-16")] <- df$datetime[df$class=="B"&df$date==as.Date("2023-01-16")] %m-% hours(1)
df %>% filter(date == as.Date("2023-01-16")) %>% ggplot(aes(x = datetime, y = CO2, color = class)) + geom_line()

to_hhmm <- function(t) {
  h <- hour(t)
  if (nchar(h) == 1) {
    h <- paste0("0",h)
  }
  m <- floor(minute(t)/10) * 10
  if (nchar(m) == 1) {
    m <- "00"
  }
  hhmm <- paste0(h, ":", m)
}

# Filter time when people are in class
rtime <- readRDS("data-clean/time-in-room.rds") 
rtime$hhmm <- sapply(rtime$datetime, to_hhmm)
 
df$hhmm <- sapply(df$datetime, to_hhmm) 
df <- df %>%
  left_join(rtime %>% select(class, date, hhmm, no_palas)) %>%
  mutate(palas = ifelse(is.na(no_palas), T, F)) 
df_filt <- df %>%
  filter(!no_palas)

saveRDS(df, "data-clean/palas-all.rds")


#### Mean particle concentration ####

df_daily_av <- df_filt %>%
  group_by(class, date) %>%
  summarize(across(c(Cn, PMtot, PM1, PM25, PM4, PM10, Temp, rH, CO2), ~ mean(.x))) %>%
  ungroup()


#### Air change rate ####

# subset data between 7 to 17 h
vent_dat <- df %>%
  group_by(class, date) %>%
  slice(1) %>%
  dplyr::select(class, date) %>%
  mutate(datetime = map(date, function(x) seq.POSIXt(as.POSIXct(paste(x, "06:50:00")), as.POSIXct(paste(x, "17:10:00")), by = "1 min"))) %>%
  unnest(datetime) %>%
  left_join(df %>% dplyr::select(class, datetime, CO2, no_palas)) %>%
  mutate(no_palas = ifelse(is.na(no_palas), T, no_palas)) %>%
  group_by(class, date) %>%
  arrange(datetime) %>%
  mutate(CO2 = zoo::na.approx(CO2, na.rm = F, maxgap = Inf)) %>%
  tidyr::fill(CO2, .direction = "downup") %>%
  ungroup() 


# add number of absences
epi_dat <- readRDS("data-clean/epidemiological.rds")
vent_dat <- vent_dat %>%
  left_join(epi_dat %>% dplyr::select(date, class, n_present)) %>%
  mutate(n_present = ifelse(no_palas, 0, n_present))

# linear increase/decrease 10min before/after room is supposed to be vacated
linear_impute <- function(x, y) {
  lag_y <- lag(y, 10)
  lag_y <- ifelse(is.na(lag_y), T, lag_y)
  lead_y <- lead(y, 10)
  lead_y <- ifelse(is.na(lead_y), T, lead_y)
  x_na <- ifelse(!y, x,
              ifelse(!lag_y, NA,
                     ifelse(!lead_y, NA, x)))
  x_imp <- zoo::na.approx(x_na, na.rm = F, maxgap = Inf)
  return(round(x_imp))
}
vent_dat <- vent_dat %>%
  group_by(class, date) %>%
  arrange(datetime) %>%
  mutate(n_present_imp = linear_impute(n_present, no_palas)) %>%
  ungroup() 

# inspect
# vent_dat %>%
#   ggplot(aes(x = datetime, y = n_present_imp)) +
#   geom_line() +
#   facet_wrap(class ~ date, scales = "free_x")

# add teacher during time in class (not outside)
vent_dat <- vent_dat %>% 
  mutate(n_present_imp = ifelse(!no_palas, n_present_imp + 1, n_present_imp))

# average CO2 emission rate
#' assume CO2 generation rate based on Persily et al. (2007): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5666301/
#' For students: assume 0.3 L/min (based on physical activity level met 1.6, age group 11 to 16, roughly the average of males and females)
#' For teachers: assume 0.36 L/min (based on physical activity level met 1.6, roughly the average of age groups 30 to 60 of males and females)
vent_dat <- vent_dat %>%
  mutate(nG = ifelse(!no_palas, (n_present_imp-1) * 0.3 + n_present_imp * 0.36, n_present_imp * 0.3))

# add volume 
vent_dat$V <- 233

# lead CO2
vent_dat <- vent_dat %>%
  rename(C = CO2) %>%
  group_by(class, date) %>%
  arrange(datetime) %>%
  mutate(C1 = lead(C)) %>%
  slice(-n()) %>%
  ungroup() %>%
  filter(hm(format(datetime, format = "%H:%M")) >= hm("07:00"),
         hm(format(datetime, format = "%H:%M")) <= hm("17:00"))

# add dt
vent_dat$dt <- 1/60

# estimate air exchange rate
#' based on Equation 16, page 7, Batterman et al. (2017): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5334699/
transient_mass_balance <- function(Q, C, nG, V, dt, Cr) {
  6 * 10^4 * nG / Q * (1 - exp(-Q / V * dt)) + (C - Cr) * exp(- Q / V * dt) + Cr
}
#' minimize residual sum of squares
min_rss <- function(data, par) {
  #' par[1] is the ventilation rate
  #' par[2] is the outdoor CO2 level
  with(data, sum( (C1 - transient_mass_balance(par[1], C, nG, V, dt, par[2])) ^ 2))
}

aer <- vent_dat %>%
  group_by(class, date) %>%
  nest() %>%
  ungroup() %>%
  #' assume ventilation rate to be between 0.01 (exact 0 not possible) and 400 and initialize with 1
  #' assume Cr to be between 375 and 425 and estimate as well
  mutate(results = lapply(data, function(X) optim(par = c(1, 400), fn = min_rss, data = X, lower = c(0.01, 375), upper = c(Inf, 425), method = "L-BFGS-B")),
         Q = sapply(results, function(x) x$par[1]),
         Cr = sapply(results, function(x) x$par[2]),
         AER = Q / 233) 

# inspect
aer %>%
  ggplot(aes(x = date, y = AER, color = class)) +
  geom_point() +
  geom_smooth()

#### Combine data ####

env_dat <- left_join(
  aer %>% dplyr::select(class, date, AER),
  df_daily_av %>% dplyr::select(class, date, Cn, PMtot, PM1, PM25, PM4, PM10, Temp, rH, CO2)
) %>%
  mutate(AER = ifelse(is.na(CO2), NA, AER),
         aircleaner = ifelse(class == "A", 
                             ifelse(date <= as.Date("2023-01-28"), "yes",
                                    ifelse(date >= as.Date("2023-02-27"), "yes", "no")),
                             ifelse(date <= as.Date("2023-01-28"), "no",
                                    ifelse(date >= as.Date("2023-02-27"), "no", "yes"))),
         weekday = weekdays(date)) %>%
  dplyr::select(class, date, aircleaner, CO2, AER, Cn, PMtot, PM1, PM25, PM4, PM10, Temp, rH)
  
saveRDS(env_dat, "data-clean/environmental.rds")


#### Check CO2 ####

# co2_time_pl <- df %>%
#   mutate(time = as.POSIXct(time, format = "%H:%M:%S")) %>%
#   ggplot(aes(x = time, y = CO2, color = class)) + 
#   geom_line() +
#   geom_vline(aes(xintercept = as.POSIXct("2023-03-14 07:00:00")), linetype = "dashed") +
#   geom_vline(aes(xintercept = as.POSIXct("2023-03-14 17:00:00")), linetype = "dashed") +
#   geom_hline(aes(yintercept = 500), linetype = "dashed", color = "red") +
#   scale_x_datetime(labels = scales::time_format(format = "%H:%M")) +
#   facet_wrap(~ date, ncol = 5) +
#   theme_bw() +
#   theme(text = element_text(size = 8), axis.title = element_blank(),
#         legend.position = "top", legend.direction = "horizontal")
# co2_time_pl
# ggsave(plot = co2_time_pl, filename = "results/checks/co2-timeseries.pdf", width = 21 / cm(1), height = 28 / cm(1))
# 
# co2_daily_pl <- df_full %>%
#   mutate(weekday = weekdays(date),
#          no_palas = ifelse(no_palas, "No School", "School")) %>%
#   filter(!(weekday %in% c("Saturday", "Sunday")),
#          hm(hhmm) >= hm("07:00"),
#          hm(hhmm) <= hm("17:00")) %>%
#   ggplot(aes(x = hhmm, y = CO2, color = no_palas)) +
#   geom_point() +
#   facet_wrap(class ~ date) 
# co2_daily_pl
# 
# df %>%
#   ggplot(aes(x = datetime, y = rH)) +
#   geom_point(shape = 21, size = .5) +
#   geom_line() +
#   labs(y = "Relative humidity (%)") +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), text = element_text(size = 8))
# 
# ggsave("results/humidity.png", width = 16 / cm(1), height = 10 / cm(1))
# 
# df_filt_2nd  <- df_filt %>%
#   select(PM1, PM25, PM4, PM10, CO2, Cn) %>%
#   set_names(c("pm1mugm", "pm25mugm", "pm4mugm", "pm10mugm", "co2ppm", "cn1cm")) %>%
#   mutate(study = "2nd")
# 
# df_filt_1st <- read_csv("../mcid/data-clean/palas.csv") %>%
#   filter(school == "School 2") %>%
#   select(c("pm1mugm", "pm25mugm", "pm4mugm", "pm10mugm", "co2ppm", "cn1cm")) %>%
#   mutate(study = "1st")
# 
# df_comp_12 <- rbind(df_filt_1st, df_filt_2nd) %>%
#   melt(c("study"))
# 
# comp_12_pl <- df_comp_12 %>%
#   ggplot(aes(color = study, y = value, x = study)) +
#   geom_boxplot() +
#   facet_wrap(~ variable, scales = "free_y") +
#   scale_color_brewer(palette = "Set1") +
#   theme(legend.position = "none", axis.title.y = element_blank()) 
# 
# comp_12_pl
# ggsave("results/comparison_studies_palas.pdf", width = 16 / cm(1), height = 12 / cm(1))
