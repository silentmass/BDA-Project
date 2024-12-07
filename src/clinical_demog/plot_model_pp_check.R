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
fit_prefix <- "faller_classification_normal-prior_0-1"

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
prior_suffix <- "normal-prior_0-1"
fits <- readRDS("models/faller_classification_normal-prior_0-1/faller_classification_models.rds")

# Plot pp_check ----

# Filter the dataframe for selected models
fit_name <- "c-PHYSICAL_SPEED_COGNITIVE_DEPRESSION"
fit <- fits[fit_name]

fits_fontsize <- 16

fits_theme <- theme(
  text = element_text(size = fits_fontsize),
  axis.title = element_text(size = fits_fontsize),
  axis.text = element_text(size = fits_fontsize),
  strip.text = element_text(size = fits_fontsize, face = "bold"),
  legend.text = element_text(size = fits_fontsize),
  legend.title = element_text(size = fits_fontsize),
  plot.background = element_rect(fill = "white", linewidth = 0),
)

title <- paste0(gsub("_", ", ", gsub("c-", " ", fit_name)), " - ", prior_suffix)

pp_plot <- pp_check(fit[[1]], ndraws = 50) +
  ggtitle(title) +
  ylim(0, 2.0) +
  theme_minimal() +
  fits_theme

filename <- paste0("all_pp_checks", "_", prior_suffix, "_", fit_name, ".png")

ggsave(paste0(c(plots_path, filename), collapse = "/"), pp_plot, width = 10, height = 7)
