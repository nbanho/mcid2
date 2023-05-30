# libraries
library(tidyverse)
library(lubridate)

# raw data
covid <- read_csv("data-raw/secondary/COVID19Test_geoRegion_w.csv") 
ili <- readxl::read_excel("data-raw/secondary/inzidenz-grippeverdacht-de.xlsx")

# preprocessing covid 
covid <- covid %>%
  filter(geoRegion == "SO") %>%
  dplyr::select(datum, anteil_pos) %>%
  mutate(year = as.numeric(substr(datum, 1, 4)),
         week = as.numeric(substr(datum, 5, 6))) %>%
  filter(year == 2023) %>%
  dplyr::select(-datum)

start_week <- as.numeric(paste0("2023", isoweek(as.Date("2023-01-01"))))
end_week <- as.numeric(paste0("2023", isoweek(as.Date("2023-03-11"))))
date <- seq.Date(as.Date("2023-01-01"), as.Date("2023-03-11"), by = "1 day")

covid_share_pos <- data.frame(date = date) %>%
  mutate(year = year(date),
         week = isoweek(date)) %>%
  left_join(covid, by = c("year", "week")) %>%
  rename(covid_prop_pos = anteil_pos)
  

# preprocessing
ili <- ili[-c(1,2),12:14] %>%
  mutate_all(as.numeric) %>%
  set_names(c("year", "week", "ili_consult")) 


# combine
commun_trans <- left_join(covid_share_pos, ili, by = c('year', 'week')) 

# save
saveRDS(commun_trans, "data-clean/secondary.rds")
