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

#otin nyt malliin parametreja, joilla näytti plotin perusteella olevan vaikutusta
fall_pooled_formula <- bf(
  YEAR_FALL ~ 1 + MMSE + TUG + EFI_EXEC_FUNC_INDEX + GCS_NEUROTRAX,
  family = "gaussian"
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

pp_check(fall_pooled_fit)


plot(data$AGE, data$FALLER)
plot(data$EFI_EXEC_FUNC_INDEX, data$FALLER)
plot(data$GCS_NEUROTRAX, data$FALLER)
plot(data$SIX_MONTHS_FALL, data$FALLER)
plot(data$YEAR_FALL, data$FALLER)
plot(data$ABC_TOTAL_PERCENT, data$YEAR_FALL)
plot(data$AGE, data$YEAR_FALL)
plot(data$EFI_EXEC_FUNC_INDEX, data$YEAR_FALL)
plot(data$GCS_NEUROTRAX, data$YEAR_FALL)
plot(data$SF36, data$YEAR_FALL)
plot(data$DGI, data$YEAR_FALL)
plot(data$BERG, data$YEAR_FALL)
plot(data$MMSE, data$YEAR_FALL)
plot(data$FSST, data$FALLER)
plot(data$TUG, data$YEAR_FALL)

