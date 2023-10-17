#### Libraries ####

library(tidyverse)
library(readxl)
library(lubridate)


#### Data ####

# read data
df <- read_xlsx("data-raw/molecular/ZF Resultate Airchecker 2023.xlsx", sheet = "Spuckproben") %>%
  rename(class = Klasse,
         date = Datum,
         virus = `Auto   Interpretation`) %>%
  dplyr::select(class, date, virus) %>% # We set A as the "Studienklasse" and B as the "Kontrollklasse" --> 3a is B and 3b is A
  mutate(class = ifelse(class == "3a", "B", "A"),
         date = as.Date(gsub("2022", "2023", as.character(date))),
         virus = ifelse(virus == "SARS-CoV-2", "CoV", virus),
         virus = ifelse(virus == "Flu B", "IFB", virus))

# line list

long_list <- df %>%
  mutate(result = ifelse(virus == "-", "Negative", "Positive")) %>%
  dplyr::select(class, date, result, virus) %>%
  arrange(class, virus) %>%
  ungroup()

saveRDS(long_list, file = "data-clean/molecular-long.rds")


#### Preprocessing ####

# counts by date 
df <- df %>%
  group_by(class, date, virus) %>%
  summarize(n = n()) %>%
  ungroup()

# add 0s
expanded_df <- expand.grid(
  class = unique(df$class), 
  date = unique(df$date),
  virus = unique(df$virus))
df <- expanded_df %>%
  left_join(df, by = c("class", "date", "virus")) %>%
  mutate(n = ifelse(is.na(n), 0, n))

# total count
df_wide <- reshape2::dcast(df, class + date ~ virus) %>%
  rename(Neg = `-`) %>%
  mutate(N = Neg + AdV + IFB + HRV + MPV + PIV + CoV)

# add additional info
df_wide <- df_wide %>%
  mutate(weekday = weekdays(date),
         aircleaner = ifelse(class == "A", 
                             ifelse(date <= as.Date("2023-01-28"), "Yes",
                                    ifelse(date >= as.Date("2023-02-27"), "Yes", "No")),
                             ifelse(date <= as.Date("2023-01-28"), "No",
                                    ifelse(date >= as.Date("2023-02-27"), "No", "Yes"))),
         n_class = ifelse(class == "A", 20, 17),
         n_class = ifelse(class == "B" & date > as.Date("2023-02-11"), 18, n_class)) %>%
  dplyr::select(class, date, weekday, aircleaner, everything())


#### Save ####

saveRDS(df_wide, file = "data-clean/molecular.rds")

