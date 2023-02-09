#### Library ####

library(tidyverse)
library(reshape2)


#### Data ####

df <- read.csv("data-raw/redcap/MCIDSchulenFollowup_DATA_2023-02-09_1356.csv") %>%
  mutate(class = ifelse(grepl("Kontroll", record_id), "A", "B"),  # A is Kontroll- and B is Studienklasse
         date = as.Date(date, format = "%d.%m.%y")) 


#### In-Room time ####

rtime <- df %>%
  select(class, date, matches("morning"), matches("afternoon")) %>%
  melt(c("class", "date")) %>%
  mutate(time = stringi::stri_extract(variable, regex = "\\d{3,4}"),
         time = ifelse(nchar(time) == 3, paste0("0", time), time),
         time = as.POSIXct(time, format = "%H%M"),
         type = ifelse(grepl("___6", variable), "break",
                       ifelse(grepl("___7", variable), "outside room",
                              ifelse(grepl("___8", variable), "partially outside room",
                                     ifelse(grepl("___9", variable), "no lesson", NA)))),
         value = ifelse(value == 1, T, F)) %>%
  filter(!is.na(date),
         !is.na(type)) %>%
  select(-variable) %>%
  dcast(class + date + time ~ type) %>%
  mutate(no_palas = ifelse(`no lesson` | `outside room` | `partially outside room` | `break`, T, F),
         time = format(time, "%H:%M:%S"),
         datetime = paste(as.character(date), as.character(time)),
         datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S")) %>%
  select(class, date, time, datetime, everything())


saveRDS(rtime, "data-clean/time-in-room.rds")
