###############################
# Analysis script: Impact of Preexisting Immunity with Host factors


#load required packages
library(tidyverse) # for data wrangling and plotting
library(lubridate) # for working with dates and time
library(here) #for data loading/saving
library(pander) # for rendering R objects into Pandoc's markdown
library(gtsummary) # for creating publication-ready summary table
library(gt) # gtsave function used to save gt table 
#webshot::install_phantomjs() # need PhantomJS to export gtsummary table

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

dat_ind <- dat_elderly %>% 
	   # keep unique participant rows for each flu season
	   distinct(id, season, age, gender, race4, bmi, obesity, prior_year_vac, dose)

# check percentage of missing values in each column by flu season
aggregate(dat_ind, by=list(dat_ind$season), FUN = function(x) { sum(is.na(x))/length(x)*100 })

#create table 1 to summarize host characteristics by dose (SD vs HD)	
summarytable <- dat_ind %>% 
	select(season, age, gender, race4, bmi, obesity, prior_year_vac, dose) %>% 
	# manually change the order in the dataset, before passing to `tbl_summary`
	mutate(gender = factor(gender, levels = c("Male", "Female"),),
				 race4  = factor(race4, levels = c("White", "Black", "Hispanic", "Other")),
				 prior_year_vac = factor(prior_year_vac, levels = c("Yes", "No", "Unknown"))) %>% 
	tbl_strata(
		strata = season,
		.tbl_fun =
			~ .x %>%
			tbl_summary(
				by = dose, missing = "no",
		    label = list(age ~ "Age",
		    						 gender ~ "Gender",
		    						 race4 ~ "Race",
		    						 bmi ~ "BMI",
		    						 obesity ~ "Obesity",
		    						 prior_year_vac ~ "Prior season vaccination"))) %>%
	bold_labels() %>% 
	# customize footnote
	as_gt() %>%
	gt::tab_source_note(gt::md("Abbreviations, SD: standard dose; HD: high dose")) 

#save table to png file for later use in manuscript
summarytable_file = here("results", "tables")
summarytable %>%
	gtsave(
		"summarytable.png", expand = 10,
		path = summarytable_file
	)

# Bivariate analysis, still in working progress
## association of variables:p values
# only show p value <0.1

cont.vars <- c("age", "bmi")
cat.vars  <- c("gender", "race4", "obesity", "prior_year_vac", "dose")
names.bivariate <- c(cont.vars, cat.vars)
	
M = diag(nrow=length(names.bivariate))
colnames(M) <- row.names(M) <- names.bivariate
for (i in 1:(nrow(M)-1)){
	for ( j in (i+1) : nrow(M)){
		bivars = c(rownames(M)[i],colnames(M)[j])
		ncat = sum(bivars %in% cat.vars)
		if(ncat ==2){
			asscn=chisq.test(table(dat_elderly[,bivars]))
			#print(asscn)
			M[i,j] <- asscn$p.value
		} else if(ncat ==1 ){
			depn = which(bivars %in% cont.vars)
			asscn= lm(dat_elderly[,bivars[depn]] ~dat_elderly[,bivars[-depn]])
			M[i,j]  <- summary(asscn)$coef[2,4]
			#print(summary(asscn))
		}else{
			asscn = lm(dat_elderly[,bivars[1]] ~dat_elderly[,bivars[2]])
			#print(summary(asscn))
			M[i,j]<- summary(asscn)$coef[2,4]
		}
	}
}
M[lower.tri(M, diag = FALSE)] <- t(M)[lower.tri(M, diag = FALSE)]
diag(M) <-0
M.long <- c()
M.long <- data.frame(col = rep(colnames(M),each=nrow(M)),
										 row =  rep(colnames(M),nrow(M)),
										 assocn = c(M))
M.long %>%
	filter(assocn<0.1) %>%
	ggplot(aes(x= factor(col, level = names.bivariate),y=factor(row, level = new.names.all2), fill = assocn)) +
	geom_tile() +
	scale_fill_gradient(low = "steelblue",
											high = "white", space = "Lab")+
	labs(x = NULL, y = NULL) +
	theme(axis.text=element_text(size=15))+   theme(axis.text.x = element_text(angle = 60, hjust = 1))


######################################
#Data visualization
######################################

# Time and dose distributions

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

# save plot 
ggsave(p1, filename = here::here("results", "figures", "days-since-vac-distribution.png"), height = 8.5, width = 11)

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

# save plot 
ggsave(p2, filename = here::here("results", "figures", "dose-distribution.png"), height = 8.5, width = 11)


# Now we should look at the time vs dose distributions as well.
p3 <- dat_elderly %>%
	tidyr::drop_na(days_before_vac) %>%
	ggplot(aes(x = days_before_vac, y = dose)) +
	geom_violin(width = 0.5, alpha = 0.5) +
	geom_jitter(alpha = 0.1, aes(col = as.factor(year))) +
	theme_bw() +
	theme(legend.position = "bottom") +
	labs(
		x = "days between vaccination and September 1st of current season",
		col = "flu season"
	)

# save plot 
ggsave(p3, filename = here::here("results", "figures", "time-dose.png"), height = 8.5, width = 11)


# Only patients >= 65 y.o. for looking at relationships w/ dose.
dat_elderly_long <- dat_long %>%
	dplyr::filter(age >= 65)

# No stratification in these ones.
dat_elderly_long %>% make_vacc_strain_plots("No-Strata", indiv = F)

# Separate by sex.
dat_elderly_long %>% make_vacc_strain_plots("Gender", indiv = F, color = gender)

# Separate by season.
dat_elderly_long %>% make_vacc_strain_plots("Season", indiv = F, color = year)

# Separate by race.
dat_elderly_long %>% make_vacc_strain_plots("Race", indiv = F, color = race4)

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