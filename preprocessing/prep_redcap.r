#### Library ####

library(tidyverse)
library(reshape2)
library(lubridate)


#### Data ####

df <- read.csv("data-raw/redcap/MCIDSchulenFollowup_DATA_2023-02-09_1356.csv") %>%
  mutate(class = ifelse(grepl("Kontroll", record_id), "B", "A"),  # A is Studien- and B is Kontrollklasse
         date = as.Date(date, format = "%Y-%m-%d")) 


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
write.csv(rtime, "data-clean/time-in-room.csv", row.names = F)


#### Absences ####

cases <- df %>%
  dplyr::select(class, initialen_sus, matches("sus_")) %>%
  filter(sus_date != "") %>%
  rename(date_symptom = sus_symp_date,
         date_absent = sus_date,
         date_back = sus_date_back,
         sex = sus_sex,
         reason = sus_absenz,
         student_id = initialen_sus) %>%
  mutate(across(c(date_symptom, date_absent, date_back), as.Date),
         sex = factor(ifelse(sex == 1, "male", "female"), levels = c("male", "female")),
         reason = factor(ifelse(is.na(reason), "unknown", ifelse(reason == 1, "sickness", "other")), levels = c("sickness", "other", "unknown")),
         respiratory_infection = `sus_symp___1` + `sus_symp___2` + `sus_symp___3` + `sus_symp___4` +
                                       `sus_symp___5` + `sus_symp___6` + `sus_symp___7` + `sus_symp___8` +
                                       `sus_symp___9` + `sus_symp___10` + `sus_symp___11` + `sus_symp___12` +
                                       `sus_symp___13`,
         respiratory_infection = ifelse(respiratory_infection >= 1, "yes", "no"),
         respiratory_infection = ifelse(sus_testing == 1 & sus_tesing_res == 1, "no", respiratory_infection), # re-label one case that tested negative
         respiratory_infection = factor(respiratory_infection, levels = c("no", "yes"))) %>% 
  dplyr::select(class, student_id, date_symptom, date_absent, date_back, sex, reason, respiratory_infection) 

# check delay between symptom onset and absence
cases <- cases %>%
  mutate(delay_absence = as.numeric(difftime(date_absent, date_symptom, units = "days"))) 
View(cases)

#' there was one case with 8 day delay between symptom onset and absence
#' in the raw report it says "since Wednesday / 25.01.2023", 
#' but the absence was on the 02.02.2023, so the last Wednesday would be 01.02.2023,
#' which would only mean 1 day of delay. This seems more likely and we code it as such.

cases$date_symptom[cases$delay_absence == 8] <- as.Date("2023-02-01")

# check time absent
cases <- cases %>%
  mutate(time_absent = as.numeric(difftime(date_back, date_absent, units = "days"))) 
View(cases)

#' the following students were absent only in the morning (and return in the afternoon)
#' nevertheless, we consider them as being absent for the whole day, 
#' thus shifting their date back by one day
cases$date_back[cases$class == "B" & cases$student_id == 7 & cases$date_back == as.Date("2023-01-19")] <- as.Date("2023-01-19") + days(1)
cases$date_back[cases$class == "A" & cases$student_id == 6 & cases$date_back == as.Date("2023-01-26")] <- as.Date("2023-01-26") + days(1)
cases$date_back[cases$class == "A" & cases$student_id == 15 & cases$date_back == as.Date("2023-02-15")] <- as.Date("2023-02-15") + days(1)
cases$date_back[cases$class == "A" & cases$student_id == 19 & cases$date_back == as.Date("2023-02-23")] <- as.Date("2023-02-23") + days(1)
cases$date_back[cases$class == "A" & cases$student_id == 12 & cases$date_back == as.Date("2023-03-02")] <- as.Date("2023-03-02") + days(1)

#' the following students return on the same date of their absence
#' but in the report the original date written with pencil,
#' which was one day later, was overwritten with pen
#' we assume that the more sensible pencil date is the correct one

cases$date_back[cases$class == "B" & cases$student_id == 5 & cases$date_back == as.Date("2023-01-16")] <- as.Date("2023-01-16") + days(1)
cases$date_back[cases$class == "B" & cases$student_id == 10 & cases$date_back == as.Date("2023-01-17")] <- as.Date("2023-01-17") + days(1)

#' there is one missing case which will be set to no respiratory infection
cases$respiratory_infection[cases$reason=="unknown"] <- "no"
#' the date of absence will be imputed assuming the median duration of absence for all non respiratory infections
med_dur_no <- median(difftime(cases$date_back[cases$respiratory_infection=="no"], 
                              cases$date_absent[cases$respiratory_infection=="no"], 
                              units = "days"), 
                     na.rm = T)
cases$date_back[cases$reason=="unknown"] <- cases$date_absent[cases$reason=="unknown"] + days(as.numeric(med_dur_no))

# save long data
cases <- cases %>% 
  dplyr::select(-delay_absence, -time_absent)
saveRDS(cases, "data-clean/redcap-long.rds")

# compute absences 
absences <- cases %>%
  mutate(date = map2(date_absent, date_back, function(x, y) if (x == y) { NA } else { seq.Date(x, y - days(1), by = "1 day") })) %>% # the if else clause is not necessary anymore after the above preprocessing steps
  unnest(date) %>%
  filter(!is.na(date)) %>%
  group_by(class, date) %>%
  summarize(n_absent = n()) %>%
  ungroup()

# compute respiratory cases
resp_cases <- cases %>%
  filter(respiratory_infection == "yes") %>%
  group_by(class, date_symptom) %>%
  summarize(new_cases = n()) %>%
  ungroup() %>%
  rename(date = date_symptom)

# full join
epi_data <- data.frame(date = rep(seq.Date(as.Date("2023-01-16"), as.Date("2023-03-11"), by = "1 day"), 2)) %>%
  mutate(class = rep(c("B", "A"), each = nrow(.) / 2),
         n_class = rep(c(17, 20), each = nrow(.) / 2),
         weekday = weekdays(date),
         weekend = ifelse(weekday %in% c("Saturday", "Sunday"), 1, 0),
         vacation = ifelse(date %in% seq.Date(as.Date("2023-02-06"), as.Date("2023-02-11"), by = "1 day"), 1, 0),
         aircleaner = ifelse(class == "A", 
                             ifelse(date <= as.Date("2023-01-28"), "Yes",
                                    ifelse(date >= as.Date("2023-02-27"), "Yes", "No")),
                             ifelse(date <= as.Date("2023-01-28"), "No",
                                    ifelse(date >= as.Date("2023-02-27"), "No", "Yes"))),
         aircleaner = ifelse(weekend + vacation > 0, "No", aircleaner)) %>%
  mutate(n_class = ifelse(class == "B" & date > as.Date("2023-02-11"), 18, n_class)) %>% # class B was 17 students before and 18 after the vacation
  left_join(absences) %>%
  left_join(resp_cases) %>%
  mutate(new_cases = ifelse(is.na(new_cases), 0, new_cases),
         n_absent = ifelse(is.na(n_absent), 0, n_absent),
         n_present = n_class - n_absent) %>%
  dplyr::select(class, n_class, date, weekday, weekend, vacation, aircleaner, n_absent, n_present, new_cases)

saveRDS(epi_data, "data-clean/epidemiological.rds")


#### Quick Look ####
cases %>% 
  group_by(class, reason, respiratory_infection) %>%
  summarize(n = n())

cases %>%
  filter(respiratory_infection == "yes") %>%
  mutate(weekday = weekdays(date_symptom)) %>%
  group_by(weekday) %>%
  summarize(n = n())

cases %>%
  filter(respiratory_infection == "yes") %>%
  mutate(dd = date_absent - date_symptom) %>%
  ggplot(aes(x = dd)) +
  geom_histogram()

cases %>%
  filter(respiratory_infection == "yes") %>%
  mutate(dd = date_back - date_absent) %>%
  ggplot(aes(x = dd)) +
  geom_histogram()
  
