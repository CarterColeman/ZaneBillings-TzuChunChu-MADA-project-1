###############################
# Exploratory data analysis script: Impact of Preexisting Immunity with Host factors


#load required packages
library(tidyverse) # for data wrangling and plotting
library(lubridate) # for working with dates and time
library(here) #for data loading/saving
library(pander) # for rendering R objects into Pandoc's markdown
library(gtsummary) # for creating publication-ready summary table
library(gt) # gtsave function used to save gt table 

# source functions for plotting
source(here::here("code","analysis_code","EDA.R"))

# load data
dat_clean <- readRDS(here::here("data","processed_data","clean_data.rds"))
dat_long  <- readRDS(here::here("data","processed_data","long_data.rds"))

# only patients >= 65 y.o.
dat_elderly <- dat_clean %>%
	dplyr::filter(age >= 65)

######################################
#Data exploration/description
######################################

# create table 1 to summarize host characteristics by dose (SD vs HD)	
summarytable <- dat_elderly %>% 
	select(season, age, gender, race2, bmi, obesity, prior_year_vac, dose) %>% 
	# manually change the order in the dataset, before passing to `tbl_summary`
	mutate(gender = factor(gender, levels = c("Male", "Female"),),
				 race2  = factor(race2, levels = c("White", "Black", "Hispanic", "Other")),
				 prior_year_vac = factor(prior_year_vac, levels = c("Yes", "No", "Unknown"))) %>% 
	tbl_strata(
		strata = season,
		.tbl_fun =
			~ .x %>%
			tbl_summary(
				by = dose, missing = "no",
		    label = list(age ~ "Age",
		    						 gender ~ "Gender",
		    						 race2 ~ "Race",
		    						 bmi ~ "BMI",
		    						 obesity ~ "Obesity",
		    						 prior_year_vac ~ "Prior season vaccination"))) %>%
	bold_labels() %>% 
	# customize footnote
	as_gt() %>%
	gt::tab_source_note(gt::md("Abbreviations, SD: standard dose; HD: high dose"))

ind_elderly <- dat_elderly %>% 
  # keep unique participant rows for each flu season
  distinct(id, season, age, gender, race2, bmi, obesity, prior_year_vac, dose)


# check percentage of missing values in each column by flu season
aggregate(ind_elderly, by=list(ind_elderly$season), FUN = function(x) { sum(is.na(x))/length(x)*100 })

# save summary table 
summarytable_location = here("results","tables", "summarytable.rds")
saveRDS(summarytable, file = summarytable_location)

######################################
#Data visualization
######################################


#### Time and dose distributions ####

# So the next thing to do is examine the distribution of times for the elderly population.

p1 <- dat_elderly %>%
	tidyr::drop_na(days_before_vac) %>%
	ggplot(aes(x = days_before_vac)) +
	geom_histogram(binwidth = 7, fill = "white", col = "black") +
	theme_bw() +
	labs(
		x = "days between vaccination and September 1st of current season") +
	facet_wrap(~season)

# Here the bin size is 7 so the histogram is effectively showing how many elderly people were vaccinated per week.
# Not sure what is going on with the 2016 flu season data so it is probably best to make sure there are not any errors here.

# Now we can quickly look at the distribution of the dose as well.

p2 <- dat_clean %>%
	tidyr::drop_na(days_before_vac) %>%
	ggplot(aes(x = days_before_vac, fill = dose)) +
	geom_histogram(
		aes(y = ..density..),
		binwidth = 7, col = "black", position = "identity", alpha = 0.5
	) +
	theme_bw() +
	labs(
		x = "days between vaccination and September 1st of current season"
	) +
	facet_wrap(~season)


# Now we should look at the time vs dose distributions as well.
p3 <- dat_elderly %>%
	tidyr::drop_na(days_before_vac) %>%
	ggplot(aes(x = days_before_vac, y = dose)) +
	geom_violin(width = 0.5, alpha = 0.5) +
	geom_jitter(alpha = 0.1, aes(col = as.factor(season))) +
	theme_bw() +
	theme(legend.position = "bottom") +
	labs(
		x = "days between vaccination and September 1st of current season",
		col = "flu season"
	)


#### Relationship between host factors and titer increase ####

# density plot
gender_p <- ggplot(dat_elderly, aes(x=titerincrease, fill=gender)) +
  geom_density(alpha=.30) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  facet_wrap(~season)

race_p <- ggplot(dat_elderly, aes(x=titerincrease, fill=race2)) +
  geom_density(alpha=.30) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  facet_wrap(~season)

obesity_p <- ggplot(dat_elderly, aes(x=titerincrease, fill=obesity)) +
  geom_density(alpha=.30) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  facet_wrap(~season)

prior_year_vac_p <- ggplot(dat_elderly, aes(x=titerincrease, fill=prior_year_vac)) +
  geom_density(alpha=.30) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  facet_wrap(~season)

dose_p <- ggplot(dat_elderly, aes(x=titerincrease, fill=dose)) +
  geom_density(alpha=.30) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  facet_wrap(~season)

# violin plot with jitter
age_plot <- ggplot(dat_elderly, aes(x=titerincrease, y=age)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(aes(colour = age), width = 0.15, height = 0.15, alpha = 0.5, size = 0.85) +
  scale_colour_gradient(low="#FFA07A",high="#7226F5") + 
  scale_y_continuous(limits = c(65,85), breaks = seq(65,85,5)) +
  coord_flip() +
  facet_wrap(~season) 


bmi_plot <- 
  dat_elderly %>% 
  filter(season > 2016) %>% 
  ggplot(aes(x=titerincrease, y=bmi)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(aes(colour = bmi), width = 0.15, height = 0.15, alpha = 0.5, size = 0.85) +
  scale_colour_gradient(low="#FFA07A",high="#7226F5") + 
  coord_flip() +
  facet_wrap(~season) 

# save plot 
ggsave(plot = p1, filename = here::here("results", "figures", "days-since-vac-distribution.png"), height = 8.5, width = 11)
ggsave(plot = p2, filename = here::here("results", "figures", "dose-distribution.png"), height = 8.5, width = 11)
ggsave(plot = p3, filename = here::here("results", "figures", "time-dose.png"), height = 8.5, width = 11)

ggsave(gender_p, filename = here::here("results", "figures", "gender_plot.png"), height = 8.5, width = 11)
ggsave(race_p, filename = here::here("results", "figures", "race_plot.png"), height = 8.5, width = 11)
ggsave(obesity_p, filename = here::here("results", "figures", "obesity_plot.png"), height = 8.5, width = 11)
ggsave(prior_year_vac_p, filename = here::here("results", "figures", "prior_year_vac_plot.png"), height = 8.5, width = 11)
ggsave(dose_p, filename = here::here("results", "figures", "dose_plot.png"), height = 8.5, width = 11)
ggsave(age_plot, filename = here::here("results", "figures", "age_plot.png"), height = 8.5, width = 11)
ggsave(bmi_plot, filename = here::here("results", "figures", "bmi_plote.png"), height = 8.5, width = 11)


  

# Only patients >= 65 y.o. for looking at relationships w/ dose.
# use long format data since it's easier to generate panel plots
dat_elderly_long <- dat_long %>%
	dplyr::filter(age >= 65)

# No stratification in these ones.
dat_elderly_long %>% make_vacc_strain_plots("No-Strata", indiv = F)

# Separate by sex.
dat_elderly_long %>% make_vacc_strain_plots("Gender", indiv = F, color = gender)

# Separate by season.
dat_elderly_long %>% make_vacc_strain_plots("Season", indiv = F, color = factor(season))

# Separate by race.
dat_elderly_long %>% make_vacc_strain_plots("Race", indiv = F, color = race2)

# Separate by dose.
dat_elderly_long %>% make_vacc_strain_plots("Dose", indiv = F, color = dose)

# Separate by age.
dat_elderly_long %>% make_vacc_strain_plots("Age", indiv = F, color = age)

# Separate by BMI
dat_elderly_long %>%
	tidyr::drop_na(bmi) %>%
	make_vacc_strain_plots("BMI", indiv = F, color = bmi)

# Separate by obesity.
dat_elderly_long %>% 
	tidyr::drop_na(obesity) %>%
	make_vacc_strain_plots("Obesity", indiv = F, color = obesity)

# Color with days since vaccination
dat_elderly_long %>%
	tidyr::drop_na(days_before_vac) %>%
	make_vacc_strain_plots("Time", indiv = F, color = days_before_vac)

# Days before vacc but only elderly
dat_elderly_long %>%
	tidyr::drop_na(days_before_vac) %>%
	make_vacc_strain_plots("Elderly-Time", indiv = F, color = days_before_vac)