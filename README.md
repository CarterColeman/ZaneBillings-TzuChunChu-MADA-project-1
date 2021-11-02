# Project Summary 
The purpose of this repository is to evaluate whether vaccine dose modify the effect of pre-exisitng immunity on immnue response to the seasonal influeza vaccine. The dataset of this project was provided to Dr. Handel research group through UGA’s CIVR-HRP site, a division of the NIH CIVICs program. The data has already been partially cleaned, but a series of data wrangling and cleaning will still be required. The final dataset includes 114 participants who were 65 years old or over (eligible for high dose Influenza vaccine) from UGA cohort during the seasons between 2016 and 2020. There are total 32 variables and 884 observations. The potential outcomes of interest are **HAI titer fold change (continuous)**, **Seroconversion (categorical)** and **Seroprotection (categorical)**. For analysis, we will first look into HAI titer fold change as the main outcome. The host characteristics we will focus on are age, sex, race, number of days from September 1st of the vaccine season to the vaccination date, prior season vaccination, and BMI.


- The project proposal is available in the `products` folder.  
- All the figures and tables were saved in the `result` folder.
- The *Manuscript.html* can be knitted from the *Manuscript.R* script located in the `products` folder.


# Description of datafiles and R codes
**Data folder:** 
- *new_clean_data.Rda* is the data downloaded directly from a private GitHub repository, and is used for the analysis
- *clean_data.Rda* and *long_data.Rda* were generated using using the *processingscript.R*

R scripts for part two of our project are located in the `code` folder, and you will see the text and codes used for data cleaning, wrangling, explorary and visulization. All the figures and tables were saved in the `result` folder.

**Code folder:** 
- *processingscript.R* runs the data cleaning and wrangling process
- *exploration.Rmd* produces tables and plottings to explore the relation between our predictors of interests and main outcome. 
- *analysisscript.R* generated the results from several different regression models and predictive modeling.   
- *EDA.R* is the function for plotting

