#### Library ####

library(tidyverse)
library(reshape2)
library(lubridate)


#### Data ####

files_A <- list.files("data-raw/palas/Palas Wochen 1-3/E3a_3 Palas_AQGuard_SN15988/aqGuard_15988", full.names = T)
files_B <- list.files("data-raw/palas/Palas Wochen 1-3/E3b_3 Palas_AQGuard_SN15973/aqGuard_15973", full.names = T)
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

df_filt <- df 
df_filt$hhmm <- sapply(df_filt$datetime, to_hhmm) 
df_filt <- df_filt %>%
  left_join(rtime %>% select(class, date, hhmm, no_palas)) %>%
  filter(!no_palas)

saveRDS(df, "data-clean/palas.rds")


#### Quick look ####

df_filt_2nd  <- df_filt %>%
  select(PM1, PM25, PM4, PM10, CO2, Cn) %>%
  set_names(c("pm1mugm", "pm25mugm", "pm4mugm", "pm10mugm", "co2ppm", "cn1cm")) %>%
  mutate(study = "2nd")

df_filt_1st <- read_csv("../mcid/data-clean/palas.csv") %>%
  filter(school == "School 2") %>%
  select(c("pm1mugm", "pm25mugm", "pm4mugm", "pm10mugm", "co2ppm", "cn1cm")) %>%
  mutate(study = "1st")

df_comp_12 <- rbind(df_filt_1st, df_filt_2nd) %>%
  melt(c("study"))

comp_12_pl <- df_comp_12 %>%
  ggplot(aes(color = study, y = value, x = study)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none", axis.title.y = element_blank()) 

comp_12_pl
ggsave("results/comparison_studies_palas.pdf", width = 16 / cm(1), height = 12 / cm(1))
