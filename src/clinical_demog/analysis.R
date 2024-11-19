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

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(patchwork)) {
  install.packages("patchwork")
  library(patchwork)
}

if(!require(cmdstanr)){
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  library(cmdstanr)
}

cmdstan_installed <- function(){
  res <- try(out <- cmdstanr::cmdstan_path(), silent = TRUE)
  !inherits(res, "try-error")
}

if(!cmdstan_installed()){
  install_cmdstan()
}

data <- readRDS("data/clinical_demog_clean.rds")
head(data)

#tässä on ensimmäinen malli, ei näytä toimivan, eli formulaa pitää muuttaa. 
#nyt siinä on vain yksi parametri (AGE), jolla ei summaryn perusteella ole vaikutusta
#periaatteessa olisi helppo tehdä kaksi eri mallia samoilla parametreilla: 
#toinen pooled ja toinen hierarchical
fall_pooled_formula <- bf(
  FALLER | trials(1) ~ 1 + AGE,
  family = binomial()
  )

fall_pooled_fit <- brm(
  formula = fall_pooled_formula,
  data = data
)

summary(fall_pooled_fit)

source("src/clinical_demog/utils/analysis_plotting_helper_functions.R")
source("src/clinical_demog/utils/analysis_helper_functions.R")

# Generate diagnostic plots
plot_mcmc_diagnostics(fall_pooled_fit)

# Check convergence metrics
check_convergence(fall_pooled_fit)