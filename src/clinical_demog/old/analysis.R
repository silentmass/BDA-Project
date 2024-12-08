# Load packages and set theme ----

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

if (!require(pROC)) {
  install.packages("pROC")
  library(pROC)
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

ggplot2::theme_set(ggplot2::theme_minimal())


# Source helper functions ----

source("src/clinical_demog/utils/analysis_plotting_helper_functions.R")
source("src/clinical_demog/utils/analysis_helper_functions.R")


# Load and preview data ----

data <- readRDS("data/clinical_demog_clean.rds")
head(data)
names(data)

# Plot data variables ----

# plot(data$AGE, data$FALLER)
# plot(data$EFI_EXEC_FUNC_INDEX, data$FALLER)
# plot(data$GCS_NEUROTRAX, data$FALLER)
# plot(data$SIX_MONTHS_FALL, data$FALLER)
# plot(data$YEAR_FALL, data$FALLER)
# plot(data$ABC_TOTAL_PERCENT, data$YEAR_FALL)
# plot(data$AGE, data$YEAR_FALL)
# plot(data$EFI_EXEC_FUNC_INDEX, data$YEAR_FALL)
# plot(data$GCS_NEUROTRAX, data$YEAR_FALL)
# plot(data$SF36, data$YEAR_FALL)
# plot(data$DGI, data$YEAR_FALL)
# plot(data$BERG, data$YEAR_FALL)
# plot(data$MMSE, data$YEAR_FALL)
# plot(data$FSST, data$FALLER)
# plot(data$TUG, data$YEAR_FALL)
# plot(data$PASE, data$YEAR_FALL)
# plot(data$S3_VELOCITY, data$YEAR_FALL)
# plot(data$BASE_VELOCITY, data$YEAR_FALL)

# Gaussian YEAR_FALL =================================

# Gaussian YEAR_FALL - MMSE + TUG + EFI_EXEC_FUNC_INDEX + GCS_NEUROTRAX ----

# otin nyt malliin parametreja, joilla näytti plotin perusteella olevan vaikutusta
fall_pooled_formula <- bf(YEAR_FALL ~ 1 +
                            MMSE +
                            TUG +
                            EFI_EXEC_FUNC_INDEX +
                            GCS_NEUROTRAX,
                          family = "gaussian")

# Show default priors
get_prior(fall_pooled_formula, data = data)

fall_pooled_fit <- brm(
  formula = fall_pooled_formula,
  data = data
)

summary(fall_pooled_fit)

# Create default density plot
pp_check(fall_pooled_fit)


# Classify subjects - FALLER - Bernoulli =================================


######## BERG + ABC + TUG ----

selected_variables <- c(
  "BERG",
  "ABC_TOTAL_PERCENT",
  "TUG"
)
fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(selected_variables, collapse = " + "))), family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(-0.5, 0.5), class = "b", coef = "BERG"),
  prior(lognormal(0, 1), class = "b", coef = "ABC_TOTAL_PERCENT"),
  prior(lognormal(0, 1), class = "b", coef = "TUG")
)

fall_class_fit <- brm(formula = fall_class_formula, data = data, family = bernoulli(), prior = fall_class_priors)

plot_pp_check(fall_class_fit, selected_variables)


######## BERG + velocities ----

selected_variables <- c(
  "BERG",
  "S3_VELOCITY",
  "BASE_VELOCITY"
)
fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(selected_variables, collapse = " + "))), family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(-0.5, 0.5), class = "b", coef = "BERG"),
  prior(normal(-0.5, 0.5), class = "b", coef = "BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "S3_VELOCITY")
)

fall_class_fit <- brm(formula = fall_class_formula, data = data, family = bernoulli(), prior = fall_class_priors)

plot_pp_check(fall_class_fit, selected_variables)


######## BERG (scaled) + velocities ----

data$BERG_std <- scale(data$BERG)[,1]

selected_variables <- c("BERG_std", "S3_VELOCITY", "BASE_VELOCITY")
fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(selected_variables, collapse = " + "))), family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(-0.5, 0.5), class = "b", coef = "BERG_std"),
  prior(normal(-0.5, 0.5), class = "b", coef = "BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "S3_VELOCITY")
)

fall_class_fit <- brm(formula = fall_class_formula, data = data, family = bernoulli(), prior = fall_class_priors)

plot_pp_check(fall_class_fit, selected_variables)


######## BERG only ----

selected_variables <- c(
  "BERG"
)
fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(selected_variables, collapse = " + "))), family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(-1, 1), class = "b", coef = "BERG")
)

fall_class_fit <- brm(formula = fall_class_formula, data = data, family = bernoulli(), prior = fall_class_priors)

plot_pp_check(fall_class_fit, selected_variables)


######## Depression: GDS + GCS_NEUROTRAX ----

selected_variables <- c(
  "GDS",
  "GCS_NEUROTRAX"
)
fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(selected_variables, collapse = " + "))), family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "GDS"),
  prior(normal(0, 1), class = "b", coef = "GCS_NEUROTRAX")
)

fall_class_fit <- brm(formula = fall_class_formula, data = data, family = bernoulli(), prior = fall_class_priors)

plot_pp_check(fall_class_fit, selected_variables)


######## Depression + BERG + velocities ----

selected_variables <- c(
  "GDS",
  "GCS_NEUROTRAX",
  "BERG",
  "BASE_VELOCITY",
  "S3_VELOCITY"
)
fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(selected_variables, collapse = " + "))), family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "GDS"),
  prior(normal(0, 1), class = "b", coef = "GCS_NEUROTRAX"),
  prior(normal(-0.5, 0.5), class = "b", coef = "BERG"),
  prior(normal(-0.5, 0.5), class = "b", coef = "BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "S3_VELOCITY")
)

fall_class_fit <- brm(formula = fall_class_formula, data = data,family = bernoulli(), prior = fall_class_priors)

plot_pp_check(fall_class_fit, selected_variables)

# All variables ----

selected_variables <- c(
  "GDS",
  "AGE",
  "GCS_NEUROTRAX",
  "EFI_EXEC_FUNC_INDEX",
  "GENDER",
  "ABC_TOTAL_PERCENT",
  "SF36",
  "MMSE",
  "MOCA",
  "FAB",
  "TUG",
  "FSST",
  "BERG",
  "DGI",
  "TMT_A",
  "TMT_B",
  "BASE_VELOCITY",
  "S3_VELOCITY",
  "PASE",
  "FEET_CLOSE_EYES_OPEN",
  "FEET_CLOSE_EYES_CLOSED",
  "TANDEM_EYES_OPEN",
  "TANDEM_EYES_CLOSED"
)
fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(selected_variables, collapse = " + "))), family = "bernoulli")

get_prior(fall_class_formula, data = data)


fall_class_priors <- c(
  # Intercept: weakly informative prior
  prior(normal(0,10), class = "Intercept"),
  
  
  # BERG: negative association with falls (higher score = better balance)
  prior(normal(0, 10), class = "b", coef = "GDS"),
  prior(normal(0, 10), class = "b", coef = "GCS_NEUROTRAX"),
  prior(normal(0, 10), class = "b", coef = "AGE"),
  prior(normal(0, 10), class = "b", coef = "EFI_EXEC_FUNC_INDEX"),
  prior(normal(0, 10), class = "b", coef = "GENDER"),
  prior(normal(0, 10), class = "b", coef = "ABC_TOTAL_PERCENT"),
  prior(normal(0, 10), class = "b", coef = "SF36"),
  prior(normal(0, 10), class = "b", coef = "MMSE"),
  prior(normal(0, 10), class = "b", coef = "MOCA"),
  prior(normal(0, 10), class = "b", coef = "FAB"),
  prior(normal(0, 10), class = "b", coef = "TUG"),
  prior(normal(0, 10), class = "b", coef = "FSST"),
  prior(normal(0, 10), class = "b", coef = "BERG"),
  prior(normal(0, 10), class = "b", coef = "DGI"),
  prior(normal(0, 10), class = "b", coef = "TMT_A"),
  prior(normal(0, 10), class = "b", coef = "TMT_B"),
  prior(normal(0, 10), class = "b", coef = "BASE_VELOCITY"),
  prior(normal(0, 10), class = "b", coef = "S3_VELOCITY"),
  prior(normal(0, 10), class = "b", coef = "PASE"),
  prior(normal(0, 10), class = "b", coef = "FEET_CLOSE_EYES_OPEN"),
  prior(normal(0, 10), class = "b", coef = "FEET_CLOSE_EYES_CLOSED"),
  prior(normal(0, 10), class = "b", coef = "TANDEM_EYES_CLOSED"),
  prior(normal(0, 10), class = "b", coef = "TANDEM_EYES_OPEN")

)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = data,
  family = bernoulli(),
  prior = fall_class_priors
)

plot_pp_check(fall_class_fit, selected_variables, variables_per_line = 5)



# All variables: Get predicted probabilities ----

summary(fall_class_fit)


posterior <- as_draws_df(fall_class_fit)
mcmc_areas(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY',
           "b_PASE","b_FEET_CLOSE_EYES_OPEN","b_FEET_CLOSE_EYES_CLOSED","b_TANDEM_EYES_OPEN","b_TANDEM_EYES_CLOSED"),  # Valitse tarkasteltavat parametrit
)


mcmc_trace(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY',
           "b_PASE","b_FEET_CLOSE_EYES_OPEN","b_FEET_CLOSE_EYES_CLOSED","b_TANDEM_EYES_OPEN","b_TANDEM_EYES_CLOSED"),  # Valitse tarkasteltavat parametrit
  prob = 0.95  # Näytä 95% uskottavuusväli
)
mcmc_intervals(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY',
           "b_PASE","b_FEET_CLOSE_EYES_OPEN","b_FEET_CLOSE_EYES_CLOSED","b_TANDEM_EYES_OPEN","b_TANDEM_EYES_CLOSED"),  # Valitse tarkasteltavat parametrit
)
## All variables with scaled variables ----

new_data <- data
new_data[, selected_variables] <- lapply(data[, selected_variables], scale)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = new_data,
  family = bernoulli(),
  prior = fall_class_priors
)

pp_check(fall_class_fit)


# All variables: Get predicted probabilities ----

summary(fall_class_fit)


posterior <- as_draws_df(fall_class_fit)
mcmc_areas(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY',
           "b_PASE","b_FEET_CLOSE_EYES_OPEN","b_FEET_CLOSE_EYES_CLOSED","b_TANDEM_EYES_OPEN","b_TANDEM_EYES_CLOSED"),  # Valitse tarkasteltavat parametrit
)


mcmc_trace(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY'),  # Valitse tarkasteltavat parametrit
  prob = 0.95  # Näytä 95% uskottavuusväli
)
mcmc_intervals(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY'),  # Valitse tarkasteltavat parametrit
  prob = 0.95  # Näytä 95% uskottavuusväli
)
library(loo)

# Lasketaan LOO-arvo brms-mallille
loo_result <- loo(fall_class_fit)

# Tulosta LOO-yhteenveto
print(loo_result)

## All scaled variables with more informative priors----
fall_class_priors <- c(
  # Intercept: weakly informative prior
  prior(normal(0,10), class = "Intercept"),
  prior(normal(0,5), class = "b")
)

new_data <- data
new_data[, selected_variables] <- lapply(data[, selected_variables], scale)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = new_data,
  family = bernoulli(),
  prior = fall_class_priors
)

pp_check(fall_class_fit)

summary(fall_class_fit)

posterior <- as_draws_df(fall_class_fit)
mcmc_areas(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY',
           "b_PASE","b_FEET_CLOSE_EYES_OPEN","b_FEET_CLOSE_EYES_CLOSED","b_TANDEM_EYES_OPEN","b_TANDEM_EYES_CLOSED"),  # Valitse tarkasteltavat parametrit
)


mcmc_trace(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY'),  # Valitse tarkasteltavat parametrit
  prob = 0.95  # Näytä 95% uskottavuusväli
)
mcmc_intervals(
  posterior,
  pars = c('b_GDS',
           'b_AGE','b_GCS_NEUROTRAX','b_EFI_EXEC_FUNC_INDEX','b_GENDER','b_ABC_TOTAL_PERCENT','b_SF36','b_MMSE','b_MOCA','b_FAB',
           'b_TUG','b_FSST','b_BERG','b_DGI','b_TMT_A','b_TMT_B','b_BASE_VELOCITY','b_S3_VELOCITY'),  # Valitse tarkasteltavat parametrit
  prob = 0.95  # Näytä 95% uskottavuusväli
)
library(loo)

# Lasketaan LOO-arvo brms-mallille
loo_result <- loo(fall_class_fit)

# Tulosta LOO-yhteenveto
print(loo_result)

