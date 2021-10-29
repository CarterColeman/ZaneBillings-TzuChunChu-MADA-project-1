---
title: "Model fitting"
author: "Tzu-Chun Chu and Zane Billings"
date: "10/28/2021"
output: html_document
---
	
# Load needed packages
library(tidymodels) # for model fitting
library(tidyverse) # data wrangling 
library(dotwhisker)  # for visualizing regression results
library(lme4) # for fitting multilevel models

# Load data
dat_clean <- readRDS(here::here("data","processed_data","clean_data.rds"))

# Check total number of participants
length(unique(clean.dat$id)) #2,286 pats

#######################################
## Model fitting: Bivariate Analysis ##
#######################################

# First, want to focus on the vaccine against H1N1, only patient >= 65 y.o., and only data after 2016 which included vaccine year and BMI data
dat_elderly <- dat_clean %>%
	dplyr::filter(age >= 65, strain_type == "H1N1", season > 2016)

# Check total number of participants again
length(unique(dat_elderly$id)) #348 pats


### Standard regression ###
# We will start by looking at the regression with only one variable at a time and ignoring the clustering
# Main outcome: titer increase (titerincrease) 
# Continuous predictors: 
# 1. Age (age)
# 2. Days to vaccination (days_before_vac)
# 3. Accumulative vaccination (cumulative_prior_vac)
# 4. BMI (bmi)
# 5. pre-vaccine titer (pretiter)

# Categorical predictors: 
# 1. Dose (SD vs. HD)
# 2. Race (Whitw vs. non-White)
# 3. Prior year vaccination (prior_year_vac)
# 4. Obesity

# Create new race variable (White vs. non-White) 
dat_elderly <- 
	dat_elderly %>% 
	mutate(white = ifelse(race2 == "White", "White", "Non-White"))

# First, let's set engine for performing linear regression
linear_fit <- 
	linear_reg() %>% 
	set_engine("lm")

## Fits a linear model to the continuous outcome using each continuous variable
lm_age <- 
	linear_fit %>% 
	fit(titerincrease ~ age, data = dat_elderly)

lm_daysbvac <- 
	linear_fit %>% 
	fit(titerincrease ~ days_before_vac, data = dat_elderly)

lm_accpvac <- 
	linear_fit %>% 
	fit(titerincrease ~ cumulative_prior_vac, data = dat_elderly)

lm_bmi <- 
	linear_fit %>% 
	fit(titerincrease ~ bmi, data = dat_elderly)

lm_pretiter <- 
	linear_fit %>% 
	fit(titerincrease ~ pretiter, data = dat_elderly)


# Summarize the results
tidy(lm_age)
tidy(lm_daysbvac)
tidy(lm_accpvac)
tidy(lm_bmi)
tidy(lm_pretiter)

# Independent t-tests were used to compare the mean of titer increase between groups
t.test(titerincrease ~ dose, data = dat_elderly, alternative = "two.sided")
t.test(titerincrease ~ white, data = dat_elderly, alternative = "two.sided")
t.test(titerincrease ~ prior_year_vac, data = dat_elderly, alternative = "two.sided")
t.test(titerincrease ~ obesity, data = dat_elderly, alternative = "two.sided")


## Fit the non-multilevel models with all predictors
lm_all <- 
	linear_fit %>% 
	fit(titerincrease ~ pretiter + dose + age + white + obesity + days_before_vac + cumulative_prior_vac + prior_year_vac, data = dat_elderly)

#Summarize the results
tidy(lm_all)

# generate a dot-and-whisker plot of linear regression results with all predictors
dotwhisker_all_pred <- tidy(lm_all) %>% 
	dotwhisker::dwplot(dot_args = list(size = 2, color = "black"),
										 whisker_args = list(color = "black"),
										 vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2))
dotwhisker_all_pred


#############################
## Linear Mixed Model Fit  ##
#############################

## Intercept only model
interceptonlymodel <- lmer(formula = titerincrease ~ 1 + (1|id), # grouping variable is subject ID
													 data    = dat_elderly) 

summary(interceptonlymodel) #get parameter estimates

# We can see under the Random effect that the residual variance on the id level is 0.4711 residual variance on the first level is 0.9437. This means that the intraclass correlation (ICC) is 0.4711/(0.9437+0.4711)=0.33.
# Under Fixed Effects the estimate of the intercept is 0.79.

## Primary predictor model, we just first add pre-vaccine titer as random factor
m1 <- lmer(formula = titerincrease ~ 1 + pretiter + (1 + pretiter|id), 
													 data    = dat_elderly) 

summary(m1)
#We can see that fixed regression slope for pre-vaccine titer was significant. However, we noticed that the error term for the slope of the variable pre-vaccine titer was very small (0.0000103). This likely means that 
#there is no slope variance of the per-vaccination titer variable between subjects, and therefore the random slope estimation can be dropped from the next model. 

## Primary predictor model, we will just add pre-vaccine titer as fixed factor
m2 <- lmer(formula = titerincrease ~ 1 + pretiter + (1 |id), 
					 data    = dat_elderly) 

summary(m2)
#pre-vaccination had a significant fixed effect, higher per-vaccine titer led to lower titer increase

## Add other predictors with fixed slopes and cross-level Interaction
# We would like to evaluate whether the differences in the relation between per-vaccine titer and titer increase in the subjects could be explained by the vaccine dose of that subjects.
m3 <- lmer(formula = titerincrease ~ 1 + pretiter + dose + pretiter:dose + age + white + obesity + days_before_vac + cumulative_prior_vac + prior_year_vac + (1|id), 
					 data    = dat_elderly) 

summary(m3)
#The cross-level interaction term between per-vaccine titer and dose was not significant 

# Remove interaction term
m4 <- lmer(formula = titerincrease ~ 1 + pretiter + dose + age + white + obesity + days_before_vac + cumulative_prior_vac + prior_year_vac + (1|id), 
					 data    = dat_elderly) 

summary(m4)
#From the results above, we can see that higher pre-vaccine titer was associaed with lower titer increase. When given high vaccine dose, participants had 0.36 increase in titer fold. Obesity, days before vaccination and prior year vaccine were 
#positively asscociated with increase in titer. 


# save result plots
ggsave(dotwhisker_all_pred, filename = here::here("results", "figures", "dotwhisker_all_pred.png"), height = 8.5, width = 11)
