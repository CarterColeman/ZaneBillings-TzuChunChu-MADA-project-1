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

theme_set(
  cowplot::theme_cowplot()
)

# load data
dat_clean <- readRDS(here::here("data","processed_data","clean_data.rds"))
dat_long  <- readRDS(here::here("data","processed_data","long_data.rds"))

# only patients >= 65 y.o.
dat_elderly <- dat_clean %>%
	dplyr::filter(
	  age >= 65,
	  # filter out only UGA data
	  startsWith(as.character(study), "UGA")
	 ) %>%
  # Drop factor levels that no longer exist
  dplyr::mutate(
    across(where(is.factor), forcats::fct_drop)
  )

dat_elderly_long <- dat_long %>%
  dplyr::filter(
    age >= 65,
    # filter out only UGA data
    startsWith(as.character(study), "UGA")
  ) %>%
  # Drop factor levels that no longer exist
  dplyr::mutate(
    across(where(is.factor), forcats::fct_drop)
  )

######################################
#Data exploration/description
######################################

# create table 1 to summarize host characteristics by dose (SD vs HD)	
summarytable <- dat_elderly %>% 
  select(season, age, gender, race2, bmi, obesity, prior_year_vac, dose,
         prevactiter) %>%
  dplyr::distinct() %>%
  tbl_strata(
    strata = season,
    .tbl_fun = ~ .x %>%
      tbl_summary(
        by = dose,
        missing = "no",
        label = list(
          age ~ "Age",
          gender ~ "Sex",
          race2 ~ "Race",
          bmi ~ "BMI",
          obesity ~ "Obesity",
          prior_year_vac ~ "Prior season vaccination",
          prevactiter ~ "Pre-vaccination HAI titer"
        ),
        type = list(
          age ~ "continuous",
          gender ~ "categorical",
          race2 ~ "categorical",
          bmi ~ "continuous",
          obesity ~ "categorical",
          prior_year_vac ~ "categorical",
          prevactiter ~ "continuous"
        )
      )
  ) %>%
  bold_labels() %>%
  modify_caption("Demographic characteristics of UGA Flu cohort by season.") %>%
  as_gt() %>%
  gt::tab_source_note(gt::md("SD: standard dose, HD: high dose"))

# table of outcome summaries

outcomestable <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  select(season, dose, seroconversion, seroprotection,
         postvactiter, titerincrease, subtype) %>%
  tbl_strata(
    strata = subtype,
    .combine_with = "tbl_stack",
    quiet = TRUE,
    .tbl_fun = ~ .x %>%
      tbl_strata(
        strata = season,
        .tbl_fun = ~ .x %>%
          tbl_summary(
            by = dose,
            missing = "no",
            label = list(
              seroconversion = "Seroconversion",
              seroprotection = "Seroprotection",
              postvactiter = "Post-vaccination HAI titer",
              titerincrease = "HAI titer fold change"
            ),
            type = list(
              seroconversion = "dichotomous",
              seroprotection = "dichotomous",
              postvactiter = "continuous",
              titerincrease = "continuous"
            )
          )
      )
  ) %>%
  bold_labels() %>%
  modify_caption("Outcome summaries for the UGAFlu cohort study.") %>%
  as_gt() %>%
  gt::tab_source_note(gt::md("SD: standard dose, HD: high dose"))

ind_elderly <- dat_elderly %>% 
  # keep unique participant rows for each flu season
  distinct(id, season, age, gender, race2, bmi, obesity, prior_year_vac, dose)


# check percentage of missing values in each column by flu season
aggregate(ind_elderly, by=list(ind_elderly$season), FUN = function(x) { sum(is.na(x))/length(x)*100 })

# save summary table 
summarytable_location = here("results","tables", "summarytable.Rds")
saveRDS(summarytable, file = summarytable_location)
saveRDS(outcomestable, file =  here("results","tables", "outcomestable.Rds"))

######################################
#Data visualization
######################################

# Univariate plots of age and BMI
covars_hist <- cowplot::plot_grid(
  dat_elderly_long %>%
    filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
    ggplot(aes(x = age)) +
    geom_histogram(color = "black", fill = "gray", binwidth = 1),
  dat_elderly_long %>%
    filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
    ggplot(aes(x = bmi)) +
    geom_histogram(color = "black", fill = "gray", binwidth = 1)
)
ggsave(
  plot = covars_hist,
  filename = here::here("results", "figures", "covars_hist.png"),
  width = 11.5, height = 8
)

# Univariate plots of titer increase and fold change
outcomes_plot <- cowplot::plot_grid(
  dat_elderly_long %>%
    filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
    ggplot(aes(x = postvactiter)) +
    geom_bar(color = "black", fill = "gray") +
    labs(x = "HAI post-vaccination titer") +
    scale_x_continuous(labels = seq(0, 10, 1), breaks = seq(0, 10, 1)) +
    scale_y_continuous(expand = expansion(c(0, 0.05), c(0,0))),
  dat_elderly_long %>%
    filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
    ggplot(aes(x = titerincrease)) +
    geom_bar(color = "black", fill = "gray") +
    labs(x = "HAI titer fold increase (log2)") +
    scale_x_continuous(labels = seq(-3, 7, 1), breaks = seq(-3, 7, 1)) +
    scale_y_continuous(expand = expansion(c(0, 0.05), c(0,0)))
)
ggsave(
  plot = outcomes_plot,
  filename = here::here("results", "figures", "outcomes_hist.png"),
  width = 11.5, height = 8
)

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
gender_p <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  ggplot(aes(x = titerincrease, fill = gender)) +
  geom_histogram(alpha=.30, position = "identity", color = "black",
                 aes(y = ..density..), binwidth = 1) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  scale_fill_manual(values = c("orange", "purple")) +
  facet_grid(subtype ~ season) +
  labs(x = "HAI titer increase", fill = "sex") +
  theme(legend.position = "bottom", legend.justification = "center")

race_p <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  mutate(race3 = factor(race2 == "White",
                        levels = c(TRUE, FALSE),
                        labels = c("White", "Other"))) %>%
  ggplot(aes(x = titerincrease, fill = race3)) +
  geom_histogram(alpha=.30, position = "identity", color = "black",
                 aes(y = ..density..), binwidth = 1) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  scale_fill_manual(values = c("orange", "purple")) +
  facet_grid(subtype ~ season) +
  labs(x = "HAI titer increase", fill = "race") +
  theme(legend.position = "bottom", legend.justification = "center")

obesity_p <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  ggplot(aes(x = titerincrease, fill = obesity)) +
  geom_histogram(alpha=.30, position = "identity", color = "black",
                 aes(y = ..density..), binwidth = 1) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  scale_fill_manual(values = c("orange", "purple")) +
  facet_grid(subtype ~ season) +
  labs(x = "HAI titer increase", fill = "obesity") +
  theme(legend.position = "bottom", legend.justification = "center")

prior_year_vac_p <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  ggplot(aes(x = titerincrease, fill = prior_year_vac)) +
  geom_histogram(alpha=.30, position = "identity", color = "black",
                 aes(y = ..density..), binwidth = 1) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  scale_fill_manual(values = c("orange", "purple")) +
  facet_grid(subtype ~ season) +
  labs(x = "HAI titer increase", fill = "prior vaccination status") +
  theme(legend.position = "bottom", legend.justification = "center")

dose_p <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  ggplot(aes(x = titerincrease, fill = dose)) +
  geom_histogram(alpha=.30, position = "identity", color = "black",
                 aes(y = ..density..), binwidth = 1) +
  scale_x_continuous(limits = c(-5,10), breaks = seq(-6,10,2)) +
  scale_fill_manual(values = c("orange", "purple")) +
  facet_grid(subtype ~ season) +
  labs(x = "HAI titer increase", fill = "dose") +
  theme(legend.position = "bottom", legend.justification = "center")

# violin plot with jitter
age_plot <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  ggplot(aes(x = age, y = titerincrease)) + 
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(subtype ~ season)


bmi_plot <- dat_elderly_long %>%
  filter(as.character(strains_fullname) == as.character(vaccine_component)) %>%
  ggplot(aes(x = bmi, y = titerincrease)) + 
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(subtype ~ season)

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
ggsave(bmi_plot, filename = here::here("results", "figures", "bmi_plot.png"), height = 8.5, width = 11)


# These plots won't be in the manuscript so I commented this out.
# # Only patients >= 65 y.o. for looking at relationships w/ dose.
# # use long format data since it's easier to generate panel plots
# dat_elderly_long <- dat_long %>%
# 	dplyr::filter(age >= 65)
#
# # No stratification in these ones.
# dat_elderly_long %>% make_vacc_strain_plots("No-Strata", indiv = F)
# 
# # Separate by sex.
# dat_elderly_long %>% make_vacc_strain_plots("Gender", indiv = F, color = gender)
# 
# # Separate by season.
# dat_elderly_long %>% make_vacc_strain_plots("Season", indiv = F, color = factor(season))
# 
# # Separate by race.
# dat_elderly_long %>% make_vacc_strain_plots("Race", indiv = F, color = race2)
# 
# # Separate by dose.
# dat_elderly_long %>% make_vacc_strain_plots("Dose", indiv = F, color = dose)
# 
# # Separate by age.
# dat_elderly_long %>% make_vacc_strain_plots("Age", indiv = F, color = age)
# 
# # Separate by BMI
# dat_elderly_long %>%
# 	tidyr::drop_na(bmi) %>%
# 	make_vacc_strain_plots("BMI", indiv = F, color = bmi)
# 
# # Separate by obesity.
# dat_elderly_long %>% 
# 	tidyr::drop_na(obesity) %>%
# 	make_vacc_strain_plots("Obesity", indiv = F, color = obesity)
# 
# # Color with days since vaccination
# dat_elderly_long %>%
# 	tidyr::drop_na(days_before_vac) %>%
# 	make_vacc_strain_plots("Time", indiv = F, color = days_before_vac)
# 
# # Days before vacc but only elderly
# dat_elderly_long %>%
# 	tidyr::drop_na(days_before_vac) %>%
# 	make_vacc_strain_plots("Elderly-Time", indiv = F, color = days_before_vac)