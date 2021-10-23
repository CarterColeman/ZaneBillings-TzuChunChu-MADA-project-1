###
# EDA: Impact of Preexisting Immunity with Host factors
# @author Zane Billings
# @since 2021-03-24
###

# Load packages
library(tidyverse)
library(lubridate)
library(here)
library(pander)

###############################################################################
# Constants
## REF_DATE controls what the "start" of each flu season is considered. Needs
## to be a string of the form "-MM-DD". Could make this more user friendly
## but I am the only one running this right now.
REF_DATE <- "-09-01"

## CI_TERMS controls which values are extracted from broom::tidy dataframes
## of linear models in order to plot estimates with their CI's.
CI_TERMS <- c("estimate", "conf.low", "conf.high")

###############################################################################
# Function for creating plots.
####

# Globally set gradient color scale to not be shades of blue
scale_colour_continuous <- function(...) {
  ggplot2::scale_colour_gradientn(colors = c("blue", "red"))
}

make_vacc_strain_plots <- function(.data, folder_name, indiv = TRUE,
                                   panel = TRUE, verbose = TRUE, ...) {
  # checks if a directory exists and creates it otherwise.
  # hardcoded to be a subdir of project Figures directory.
  dir.create(file.path(here::here("results","figures"), folder_name))

  # Loop through each subtype and make a panel plot for each one
  for (st in levels(.data$subtype)) {
    plt_data <- .data %>%
      dplyr::filter(
        subtype == st,
        strain_type == st
      ) %>%
      dplyr::mutate(
        across(where(is.factor), forcats::fct_drop)
      )

    # Make a sub-subdir for this subtype.
    dir.create(file.path(here::here("results","figures"), folder_name, st))

    # For each vaccine component, make a panel plot
    for (vac in levels(plt_data$vaccine_component)) {

      # Make individual plots
      if (isTRUE(indiv)) {
        for (strain in levels(plt_data$strains_fullname)) {
          # Plot log2increase vs log2pretiter
          p <- plt_data %>%
            dplyr::filter(
              vaccine_component == vac,
              strains_fullname == strain
            ) %>%
            ggplot(aes(x = prevactiter, y = titerincrease, ...)) +
            geom_jitter(width = 0.15, height = 0.15) +
            geom_smooth(method = "lm", formula = "y ~ x") +
            ggtitle(paste0("Vaccine: ", vac, "\n", "Strain: ", strain)) +
            ggpubr::stat_regline_equation(label.x = 3) +
            cowplot::theme_cowplot()
          # Color homologous plots to make easier to tell apart
          if (vac == strain) {
            p <- p + theme(panel.background = element_rect(fill = "papayawhip"))
          }
          # Save plot
          ggsave(
            filename = here::here(
              "results",
              "figures",
              folder_name,
              st,
              paste0(vac, "_", strain, ".png")
            ),
            plot = p,
            width = 7,
            height = 7,
            units = "in"
          )
          if (isTRUE(verbose)) {
            cat("saved indiv plot\n")
          }
        } # end of strain loop
      } # end of indiv plot code
      
      # Make panel plots
      if (isTRUE(panel)) {
        panel_plot <- plt_data %>%
          dplyr::filter(vaccine_component == vac) %>%
          dplyr::mutate(
            strain_name = ifelse(
              vac != "None",
              forcats::fct_relevel(strains_fullname, as.character(vac)),
              strains_fullname
            )
          ) %>%
          ggplot(aes(x = prevactiter, y = titerincrease, ...)) +
          geom_jitter(width = 0.15, height = 0.15, alpha = 0.5, size = 0.85) +
          geom_smooth(method = "lm", formula = "y ~ x", se = TRUE) +
          ggpubr::stat_regline_equation(label.x = 3) +
          ggtitle(paste0("Vaccine: ", vac)) +
          facet_wrap(~strains_fullname) +
          cowplot::theme_cowplot()
        
        ggsave(
          filename = here::here(
            "results",
            "figures",
            folder_name,
            st,
            paste0("panel_", vac, ".png")
          ),
          plot = panel_plot,
          width = 12,
          height = 7,
          units = "in"
        )
        if (isTRUE(verbose)) {
          cat("saved panel plot\n")
        }
      } # end of panel plot code
    } # end of vac loop
  } # end of st loop
} # end of function

###############################################################################
# Summary coeff data
###

make_model <- function(.data, covariates = NULL) {
  model_formula <- paste0("titerincrease ~ prevactiter")
  if (!is.null(covariates)) {
    covariates <- paste0(covariates, collapse = " + ")
    model_formula <- paste(model_formula, covariates, sep = " + ")
  }
  model <- lm(formula = as.formula(model_formula),
              data = .data)
  return(model)
}



# TODO alter this function so that it accepts arguments for strata.
test_plot_rename_this <- function() {
  dat_models %>%
    dplyr::filter(
      vaccine_component != "None",
      isFALSE(homologous),
      dose == "SD"
    ) %>%
    ggplot(aes(x = year_diff, y = estimate, col = vaccine_component)) +
    geom_point() +
    geom_smooth(method = "lm", formula = "y ~ x") +
    labs(
      x = "# years between test strain and vaccine strain circulation",
      y = "prevactiter slope coefficient"
    ) +
    facet_wrap(~subtype) +
    theme_bw()
}

# TODO adapt this code to generate this plot for different strata?
# TODO fix this code to include pre-B strains.

# Pointrange plots that summarize linear regression outputs
# TODO update this to accept different strata.
# TODO update this to save each plot to file rather than pasting together.
create_pointrange_coef_plot <- function(.data) {
  plot_list <- list()
  for (st in unique(.data$subtype)) {
    p <- .data %>%
      dplyr::filter(subtype == st, vaccine_component != "None") %>%
      ggplot(aes(x = estimate, y = strains_fullname, col = homologous)) +
      geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
      geom_vline(xintercept = 0, lty = 2, col = "firebrick3") +
      facet_wrap(~vaccine_component) +
      theme_bw() +
      labs(
        x = "prevactiter slope coefficient", y = "", title = paste(st)
      ) +
      theme(legend.position = "none")
    plot_list <- append(plot_list, list(p))
  }
  cowplot::plot_grid(plotlist = plot_list)
}

 ###############################################################################
# Function to remake all wanted individual-level scatterplots
# TODO modify this to pass indiv and panel arguments to all constituents
# TODO potentially modify this to control which plots get made? 
#  but that may be too complicated vs just running the plots.
make_all_indiv_and_panel_plots <- function() {
###
# Make the figures similar to Figure 3 for H1N1 vaccines in the reference.
###
  
  # No stratification in these ones.
  dat_long %>% make_vacc_strain_plots("No-Strata")
  
  # Separate by sex.
  dat_SD_only %>% make_vacc_strain_plots("Gender", color = gender)
  
  # Separate by season.
  dat_SD_only %>% make_vacc_strain_plots("Season", color = year)
  
  # Separate by race.
  dat_SD_only %>% make_vacc_strain_plots("Race", color = race2)
  
  # Separate by dose.
  dat_elderly %>% make_vacc_strain_plots("Dose", color = dose)
  
  # Separate by age.
  dat_SD_only %>% make_vacc_strain_plots("Age", color = age)
  
  # Separate by BMI--UGA data only
  dat_SD_only %>%
    tidyr::drop_na(bmi) %>%
    make_vacc_strain_plots("BMI", color = bmi)
  
  # Color with days since vaccination
  dat_SD_only %>%
    tidyr::drop_na(days_before_vac) %>%
    make_vacc_strain_plots("Time", color = days_before_vac)
  
  # Days before vacc but only elderly
  dat_elderly %>%
    tidyr::drop_na(days_before_vac) %>%
    make_vacc_strain_plots("Elderly-Time", color = days_before_vac)

###
# Plots with two stratifying variables.
###
  
  # Dose and age
  dat_elderly %>%
    make_vacc_strain_plots("Age-Dose", color = age, shape = dose, lty = dose)
  
  # Dose and gender
  dat_elderly %>%
    make_vacc_strain_plots("Gender-Dose", color = gender, shape = dose,
                           lty = dose)
  
  # Gender and age
  dat_SD_only %>%
    make_vacc_strain_plots("Age-Gender", color = age, shape = gender,
                           lty = gender)
  
  # Age and BMI??
  dat_SD_only %>%
    tidyr::drop_na(bmi) %>%
    make_vacc_strain_plots("Age-BMI", color = age, size = bmi)
  
  # Dose and year
  dat_elderly %>%
    make_vacc_strain_plots("Season-Dose", color = year, shape = dose,
                           lty = dose)
  
  # Dose and race
  dat_elderly %>%
    # Not enough observations in each group to look at Hispanic and Other 
    # categories
    dplyr::mutate(
      race3 = forcats::fct_other(
        race2, keep = "white", other_level = "not white"
      )
    ) %>%
    make_vacc_strain_plots("Race-Dose", color = race2, shape = dose,
                           lty = dose)
  
  # Days before vacc AND dose but only elderly
  dat_elderly %>%
    tidyr::drop_na(days_before_vac) %>%
    make_vacc_strain_plots("Elderly-Time-Dose", 
                           col = days_before_vac, shape = dose, lty = dose)
  
} # end of plotting function


