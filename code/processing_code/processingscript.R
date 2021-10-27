###############################
# processing script
###############################

#this script loads the cleaned data from the CIVIC repo, cleans it with additional steps
#and saves it as Rds file in the processed_data folder

#load required packages
library(here) #for data loading/saving
library(tidyverse) # for data wrangling
library(lubridate) # working with date and time data 

# load data
dat_orig <- readRDS(here::here("data","raw_data","new_clean_data.rds"))


######################################
# Data cleaning and wrangling
######################################

# look up data structure and summary
str(dat_orig)
summary(dat_orig)

# check percentage of missing values in each column by flu season
aggregate(dat_orig[,3:30], by=list(dat_orig$season), FUN = function(x) { sum(is.na(x))/length(x)*100 })

# It looks like 'bmi' and 'date_vaccinated' was not collected in 2014, 2015, and were mostly missing in 2016. Large proportion of 'race2'   
# was missing in 2019 and some were missing in 2020. Need to recode 'race' column. 

# Fixing an issue with race2 variable--this was fixed upstream.
# # check race variables
# table(dat_orig$race)
# table(dat_orig$race2)
# table(dat_orig$race, dat_orig$race2)


# white <- c("Polynesian","White","White/ Caucasian","White/Caucaian","White/Caucasian")
# black <- c("Black","Black or African American","Black/African American")
# hisp  <- unique(dat_orig$race)[grep("Hispanic", unique(dat_orig$race))] # extract values that contain "Hispanic"

dat_clean <- dat_orig %>% 
  # Create a categorical variable for obesity.
	mutate(obesity = ifelse(bmi >= 30, "Yes", "No")) %>%
  # Rename "days_since_vac" for backwards compatibility with analysis code
  # (due to updating dataset after code was written)
  rename(days_before_vac = days_since_vac)
  
# create days before vaccination variable -- this was also done upstream,
# after this code was already written.
# REF_DATE <- "-09-01"

# dat_clean2 <- dat_clean %>%
# 	tidyr::drop_na(prevactiter, titerincrease) %>%
# 	dplyr::mutate(
# 		year = factor(season),
# 		days_before_vac = lubridate::time_length(
# 			x = ifelse(
# 				test = is.na(season),
# 				yes = NA,
# 				no = lubridate::interval(
# 					start = lubridate::ymd(paste0(season, REF_DATE)),
# 					end = date_vaccinated
# 				)
# 			),
# 			unit = "day"
# 		)
# 	)

# pivot the data. This makes it easier to plot for all strains.
dat_long <- dat_clean %>%
	tidyr::pivot_longer(
		ends_with("vaccine_fullname"),
		names_to = "subtype", values_to = "vaccine_component"
	) %>%
	dplyr::mutate(
		subtype = factor(case_when(
			startsWith(subtype, "h1n1") ~ "H1N1",
			startsWith(subtype, "h3n2") ~ "H3N2",
			startsWith(subtype, "bvic") ~ "B-Vic",
			startsWith(subtype, "yama") ~ "B-Yam",
			TRUE ~ "uh-oh"
		))
	)

# Save as RDS
saveRDS(dat_clean, here::here("data","processed_data","clean_data.rds"))
saveRDS(dat_long, here::here("data","processed_data", "long_data.Rds"))
