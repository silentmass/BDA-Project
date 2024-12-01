######## Run full ----

# Init ----

#### Load packages and set theme ----

if (!require(tidyverse)) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require(tidybayes)) {
  install.packages("tidybayes")
  library(tidybayes)
}

if (!require(brms)) {
  install.packages("brms")
  library(brms)
}


if (!require(bayesplot)) {
  install.packages("bayesplot")
  library(bayesplot)
}

if (!require(loo)) {
  install.packages("loo")
  library(loo)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(moments)) {
  install.packages("moments")
  library(moments)
}


if (!require(pROC)) {
  install.packages("pROC")
  library(pROC)
}

if (!require(patchwork)) {
  install.packages("patchwork")
  library(patchwork)
}

if (!require(igraph)) {
  install.packages("igraph")
  library(igraph)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require(tidyr)) {
  install.packages("tidyr")
  library(tidyr)
}

if(!require(cmdstanr)){
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  library(cmdstanr)
}

if (!require(projpred)) {
  install.packages("projpred")
  library(projpred)
}

cmdstan_installed <- function(){
  res <- try(out <- cmdstanr::cmdstan_path(), silent = TRUE)
  !inherits(res, "try-error")
}

if(!cmdstan_installed()){
  install_cmdstan()
}

ggplot2::theme_set(ggplot2::theme_minimal())

#### Set SEED ----
fit_prefix <- "faller_classification_top2"

results_path <- paste0(c("results", fit_prefix), collapse = "/")
models_path <- paste0(c("models", fit_prefix), collapse = "/")
plots_path <- paste0(c("plots", fit_prefix), collapse = "/")
mcmc_path <- paste0(c("plots", fit_prefix, "MCMC"), collapse = "/")
mcmc_areas_path <- paste0(c("plots", fit_prefix, "MCMC", "areas"), collapse = "/")
mcmc_traces_path <- paste0(c("plots", fit_prefix, "MCMC", "traces"), collapse = "/")
mcmc_intervals_path <- paste0(c("plots", fit_prefix, "MCMC", "intervals"), collapse = "/")

SEED = 2024

set.seed(SEED)


#### Source helper functions ----
utils_path <- paste0(c("src", "clinical_demog", "utils"), collapse = "/")

source(paste0(c(utils_path, "analysis_plotting_helper_functions.R"), collapse = "/"))
source(paste0(c(utils_path, "analysis_helper_functions.R"), collapse = "/"))

# Load models
fits <- readRDS("models/faller_classification_top2/faller_classification_models.rds")

# PLot MCMC ----

pars <- list(
  speed = c(
    "b_Intercept", 
    "b_z_BASE_VELOCITY"
  ),  # Note the escaped bracket
  speed_dep = c(
    "b_Intercept", 
    "b_z_BASE_VELOCITY", 
    "b_z_GDS"
  ),
  speed_hier_age_group = c(
    "b_Intercept", 
    "b_z_BASE_VELOCITY", 
    "r_AGE_GROUP[1,Intercept]",
    "r_AGE_GROUP[2,Intercept]",
    "r_AGE_GROUP[3,Intercept]",
    "r_AGE_GROUP[1,z_BASE_VELOCITY]",
    "r_AGE_GROUP[2,z_BASE_VELOCITY]",
    "r_AGE_GROUP[3,z_BASE_VELOCITY]"
  ),
  speed_dep_hier_gender = c(
    "b_Intercept", 
    "b_z_BASE_VELOCITY",
    "b_z_GDS",
    "r_GENDER[0,Intercept]",
    "r_GENDER[1,Intercept]",
    "r_GENDER[0,z_BASE_VELOCITY]",
    "r_GENDER[1,z_BASE_VELOCITY]",
    "r_GENDER[0,z_GDS]",
    "r_GENDER[1,z_GDS]"
  ),
  speed_hier_age_group_v2 = c(
    "b_Intercept", 
    "b_z_BASE_VELOCITY", 
    "r_AGE_GROUP[1,Intercept]",
    "r_AGE_GROUP[2,Intercept]",
    "r_AGE_GROUP[3,Intercept]",
    "r_AGE_GROUP[1,z_BASE_VELOCITY]",
    "r_AGE_GROUP[2,z_BASE_VELOCITY]",
    "r_AGE_GROUP[3,z_BASE_VELOCITY]"
  ),
  speed_dep_hier_gender_v2 = c(
    "b_Intercept", 
    "b_z_BASE_VELOCITY",
    "b_z_GDS",
    "r_GENDER[0,Intercept]",
    "r_GENDER[1,Intercept]",
    "r_GENDER[0,z_BASE_VELOCITY]",
    "r_GENDER[1,z_BASE_VELOCITY]",
    "r_GENDER[0,z_GDS]",
    "r_GENDER[1,z_GDS]"
  )
)

fontsize = 14
base_theme <- theme_minimal(base_size = fontsize) +
  theme(
    plot.title = element_text(hjust = 0.5, size = fontsize, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = fontsize, color = "gray40"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = fontsize),
    legend.text = element_text(size = fontsize),
    legend.title = element_text(size = fontsize),
    plot.background = element_rect(fill = "white", linewidth = 0)
  )

# Loop through each fit and create plot
plots_list <- lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  mcmc_area_plot <- mcmc_areas(fit, pars = pars[[model_name]]) +
    labs(
      title = model_name,
      x = "Standardized Coefficient",
      y = "Parameter",
      subtitle = "Shaded areas represent 95% credible intervals"
    ) +
    base_theme
  mcmc_traces_plot <- mcmc_trace(fit, pars = pars[[model_name]]) +
    labs(
      title = "MCMC Chain Traces",
      x = "Iteration",
      y = "Parameter Value",
      subtitle = model_name
    ) +
    base_theme +
    theme(
      panel.grid.major = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    facet_wrap(~parameter, ncol = 4)
  mcmc_intervals_plot <- mcmc_intervals(fit, pars = pars[[model_name]]) +
    labs(
      title = model_name,
      subtitle = "95% Credible Intervals",
      x = "Standardized Effect",
      y = "Parameter"
    ) +
    base_theme
  
  # mcmc_area_plot <- plot_mcmc(fit, model_name, all_predictors[[model_name]], "area")
  # mcmc_traces_plot <- plot_mcmc(fit, model_name, all_predictors[[model_name]], "trace")
  # mcmc_intervals_plot <- plot_mcmc(fit, model_name, all_predictors[[model_name]], "interval")
  
  ggsave(paste0(c(mcmc_path, "areas", paste0("areas_", model_name, ".png")), collapse = "/"), mcmc_area_plot, width = 20, height = 20)
  ggsave(paste0(c(mcmc_path, "traces", paste0("traces_", model_name, ".png")), collapse = "/"), mcmc_traces_plot, width = 20, height = 20)
  ggsave(paste0(c(mcmc_path, "intervals", paste0("intervals_", model_name, ".png")), collapse = "/"), mcmc_intervals_plot, width = 20, height = 20)
  return(list(name=model_name, area=mcmc_area_plot, traces=mcmc_traces_plot, intervals=mcmc_intervals_plot))
})