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
SEED = 2024

set.seed(SEED)

default_prior_b <- prior(student_t(3, 0, 2.5), class = "b")

fit_prefix <- "faller_classification_top2"

results_path <- paste0(c("results", fit_prefix), collapse = "/")
models_path <- paste0(c("models", fit_prefix), collapse = "/")
plots_path <- paste0(c("plots", fit_prefix), collapse = "/")
mcmc_path <- paste0(c("plots", fit_prefix, "MCMC"), collapse = "/")
mcmc_areas_path <- paste0(c("plots", fit_prefix, "MCMC", "areas"), collapse = "/")
mcmc_traces_path <- paste0(c("plots", fit_prefix, "MCMC", "traces"), collapse = "/")
mcmc_intervals_path <- paste0(c("plots", fit_prefix, "MCMC", "intervals"), collapse = "/")

dir.create(results_path, recursive = TRUE, showWarnings = FALSE)
dir.create(models_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)
dir.create(mcmc_path, recursive = TRUE, showWarnings = FALSE)
dir.create(mcmc_areas_path, recursive = TRUE, showWarnings = FALSE)
dir.create(mcmc_traces_path, recursive = TRUE, showWarnings = FALSE)
dir.create(mcmc_intervals_path, recursive = TRUE, showWarnings = FALSE)

#### Source helper functions ----

utils_path <- paste0(c("src", "clinical_demog", "utils"), collapse = "/")

source(paste0(c(utils_path, "analysis_plotting_helper_functions.R"), collapse = "/"))
source(paste0(c(utils_path, "analysis_helper_functions.R"), collapse = "/"))

#### Load and preview data ----

data <- readRDS("data/clinical_demog_clean.rds")
head(data)
names(data)

# data$FALLER <- factor(data$FALLER)

data <- data %>%
  mutate(AGE_GROUP = case_when(
    AGE < 75 ~ 1,
    AGE >= 75 & AGE <= 80 ~ 2,
    AGE > 80 ~ 3
  ))

# Log transform TMTs
data$log_TMT_A <- log(data$TMT_A)
data$log_TMT_B <- log(data$TMT_B)

#### Log transform TMT and z-score scale and fit ----

original_cols_list <- list(
  PHYSICAL = c("AGE", "GENDER", "DGI", "TUG", "FSST"),
  SPEED = c("BASE_VELOCITY"),
  COGNITIVE = c("GCS_NEUROTRAX", "log_TMT_B"),
  DEPRESSION = c("GDS")
)
original_predictor_categories <- get_predictor_categories_from_cols_list(original_cols_list)

# Scale all predictors including log-transformed TMT variables
data[paste0("z_", names(original_predictor_categories))] <- scale(data[names(original_predictor_categories)]) # Scale also GENDER

# Test training ----
# Weakly informative directional priors
priors_weak_dir <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(-0.5, 1), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(0, 1), class = "b")
)

# 1. Simple speed model
fit_speed <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY,
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_weak_dir,
  cores = 4
)

# 2. Physical measures + speed
fit_phys_speed <- brm(
  FALLER ~ 1 + z_AGE + z_GENDER + z_DGI + z_TUG + z_FSST + z_BASE_VELOCITY,
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_weak_dir,
  cores = 4
)

# 3. Hierarchical model with additional priors for group-level parameters
priors_hier <- c(
  priors_weak_dir,
  prior(normal(0, 0.5), class = "sd"),  # group-level SDs
  prior(lkj(2), class = "cor")          # group-level correlations
)

fit_hier_age_speed <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_hier,
  cores = 4
)

summary(fit_speed)
summary(fit_phys_speed)
summary(fit_hier_age_speed)

# Set 3 ----
# Scale-appropriate priors
priors_scale <- c(
  prior(normal(0, 0.7), class = "Intercept"),  # probit scale typically < 2
  prior(normal(-0.3, 0.7), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(0, 0.7), class = "b")  # for other predictors
)

# 1. Simple speed model
fit_speed_scale <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY,
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_scale,
  cores = 4
)

# 2. Physical measures + speed
fit_phys_speed_scale <- brm(
  FALLER ~ 1 + z_AGE + z_GENDER + z_DGI + z_TUG + z_FSST + z_BASE_VELOCITY,
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_scale,
  cores = 4
)

# 3. Hierarchical model
priors_hier_scale <- c(
  priors_scale,
  prior(normal(0, 0.3), class = "sd"),  # tighter prior for group-level SDs
  prior(lkj(2), class = "cor")
)

fit_hier_age_speed_scale <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_hier_scale,
  cores = 4
)

summary(fit_speed_scale)
summary(fit_phys_speed_scale)
summary(fit_hier_age_speed_scale)

# Gender with set 3 ----

fit_hier_gender_speed_scale <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_hier_scale,
  cores = 4
)

fit_hier_gender_speed_depression_scale <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_hier_scale,
  cores = 4
)

fit_hier_gender_speed_cog_scale <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GCS_NEUROTRAX + z_log_TMT_B + 
    (1 + z_BASE_VELOCITY + z_GCS_NEUROTRAX + z_log_TMT_B | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_hier_scale,
  cores = 4
)

summary(fit_hier_gender_speed_scale)
summary(fit_hier_gender_speed_depression_scale)
summary(fit_hier_gender_speed_cog_scale)

# Gender ----

fit_hier_age_speed_cog_scale <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GCS_NEUROTRAX + z_log_TMT_B + 
    (1 + z_BASE_VELOCITY + z_GCS_NEUROTRAX + z_log_TMT_B | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_hier_scale,
  cores = 4
)

# Model with both age and gender hierarchical effects
fit_hier_age_gender_speed_scale <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + 
    (1 + z_BASE_VELOCITY | AGE_GROUP) +
    (1 + z_BASE_VELOCITY | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = priors_hier_scale,
  cores = 4
)

summary(fit_hier_age_speed_cog_scale)
summary(fit_hier_age_gender_speed_scale)


# ----
prior_summary(fit_speed)
prior_summary(fit_phys_speed)
prior_summary(fit_hier_age_speed)
prior_summary(fit_speed_scale)
prior_summary(fit_phys_speed_scale)
prior_summary(fit_hier_age_speed_scale)
prior_summary(fit_hier_gender_speed_scale)
prior_summary(fit_hier_gender_speed_depression_scale)
prior_summary(fit_hier_gender_speed_cog_scale)
prior_summary(fit_hier_age_speed_cog_scale)
prior_summary(fit_hier_age_gender_speed_scale)

summary(fit_speed)
summary(fit_phys_speed)
summary(fit_hier_age_speed)
summary(fit_speed_scale)
summary(fit_phys_speed_scale)
summary(fit_hier_age_speed_scale)
summary(fit_hier_gender_speed_scale)
summary(fit_hier_gender_speed_depression_scale)
summary(fit_hier_gender_speed_cog_scale)
summary(fit_hier_age_speed_cog_scale)
summary(fit_hier_age_gender_speed_scale)

# Compare models using LOO
loo_compare <- loo_compare(
  loo(fit_speed), 
  loo(fit_phys_speed),
  loo(fit_hier_age_speed),
  loo(fit_speed_scale),
  loo(fit_phys_speed_scale),
  loo(fit_hier_age_speed_scale),
  loo(fit_hier_gender_speed_scale),
  loo(fit_hier_gender_speed_depression_scale),
  loo(fit_hier_gender_speed_cog_scale),
  loo(fit_hier_age_speed_cog_scale),
  loo(fit_hier_age_gender_speed_scale)
)

print(loo_compare)

# fine tune best two model priors ----

# Best model - Gender, Speed, Depression with suggested priors
fit_hier_gender_speed_depression_scale_new <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # Fixed effects
    prior(normal(0, 1), class = "Intercept"),
    prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
    prior(normal(0, 0.5), class = "b", coef = "z_GDS"),
    # Group-level effects
    prior(normal(0, 0.4), class = "sd", lb = 0),
    prior(lkj_corr_cholesky(2), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Second best model - Age, Speed with suggested priors
fit_hier_age_speed_scale_new <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # Fixed effects
    prior(normal(0, 1), class = "Intercept"),
    prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
    # Group-level effects
    prior(normal(0, 0.4), class = "sd", lb = 0),
    prior(lkj_corr_cholesky(2), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# After fitting, compare with original models
loo_compare(
  loo(fit_hier_gender_speed_depression_scale),
  loo(fit_hier_gender_speed_depression_scale_new),
  loo(fit_hier_age_speed_scale),
  loo(fit_hier_age_speed_scale_new)
)

summary(fit_hier_gender_speed_depression_scale_new)
summary(fit_hier_age_speed_scale_new)

# fine tune more ----

# Version 1: More informative priors based on domain knowledge
fit_hier_gender_speed_depression_v1 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # More informative fixed effects
    prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"), # Stronger neg prior
    prior(normal(0.3, 0.3), class = "b", coef = "z_GDS"), # Slightly positive prior
    # Tighter group-level effects
    prior(exponential(3), class = "sd"),  # More concentrated near 0
    prior(lkj_corr_cholesky(4), class = "L")  # Stronger prior toward independence
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Version 2: Wider, less informative priors
fit_hier_gender_speed_depression_v2 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # Wider fixed effects
    prior(normal(0, 1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b", coef = "z_BASE_VELOCITY"),
    prior(normal(0, 1), class = "b", coef = "z_GDS"),
    # Wider group-level effects
    prior(student_t(3, 0, 0.5), class = "sd"),  # Heavier tails
    prior(lkj_corr_cholesky(1), class = "L")  # Uniform on correlations
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Version 3: Student-t priors for robustness
fit_hier_gender_speed_depression_v3 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # Student-t fixed effects
    prior(student_t(3, 0, 0.7), class = "Intercept"),
    prior(student_t(3, -0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
    prior(student_t(3, 0, 0.5), class = "b", coef = "z_GDS"),
    # Half-Cauchy for SDs
    prior(cauchy(0, 0.3), class = "sd", lb = 0),
    prior(lkj_corr_cholesky(2), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Compare all versions
loo_compare(
  loo(fit_hier_gender_speed_depression_scale),
  loo(fit_hier_gender_speed_depression_v1),
  loo(fit_hier_gender_speed_depression_v2),
  loo(fit_hier_gender_speed_depression_v3)
)

summary(fit_hier_gender_speed_depression_scale)
summary(fit_hier_gender_speed_depression_v1)
summary(fit_hier_gender_speed_depression_v2)
summary(fit_hier_gender_speed_depression_v3)

# Recommended model 1 ----

fit_recommended <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"),
    prior(normal(0.3, 0.3), class = "b", coef = "z_GDS"),
    prior(exponential(3), class = "sd"),
    prior(lkj_corr_cholesky(4), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

loo_compare(loo(fit_hier_age_speed_scale), loo(fit_recommended))

summary(fit_hier_age_speed_scale)
summary(fit_recommended)

pairs(fit_recommended)

summary(fit_recommended)

# fine tune seconde best model ----
# Version 1: More informative priors similar to previous best model
fit_hier_age_speed_v1 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"),
    prior(exponential(3), class = "sd"),
    prior(lkj_corr_cholesky(4), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Version 2: Student-t priors for robustness
fit_hier_age_speed_v2 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(student_t(3, 0, 0.7), class = "Intercept"),
    prior(student_t(3, -0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
    prior(student_t(3, 0, 0.3), class = "sd", lb = 0),
    prior(lkj_corr_cholesky(2), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Version 3: Half-Cauchy for variance components
fit_hier_age_speed_v3 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(normal(0, 0.7), class = "Intercept"),
    prior(normal(-0.5, 0.4), class = "b", coef = "z_BASE_VELOCITY"),
    prior(cauchy(0, 0.2), class = "sd", lb = 0),
    prior(lkj_corr_cholesky(2), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Compare all versions
loo_compare(
  loo(fit_hier_age_speed_scale),
  loo(fit_hier_age_speed_v1),
  loo(fit_hier_age_speed_v2),
  loo(fit_hier_age_speed_v3)
)

summary(fit_hier_age_speed_scale)
summary(fit_hier_age_speed_v1)
summary(fit_hier_age_speed_v2)
summary(fit_hier_age_speed_v3)

# Recommended model 2 ----


fit_hier_age_speed_recommended <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"),
    prior(exponential(3), class = "sd"),
    prior(lkj_corr_cholesky(4), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.9)  # Increased from default 0.8
)

# To examine the pairs plot:
pairs(fit_hier_age_speed_recommended)

# v2 ----

fit_hier_age_speed_recommended_v2 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # Less informative prior for intercept
    prior(normal(0, 1), class = "Intercept"),
    # Keep informative prior for base velocity
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"),
    # Try half-Student-t for SD parameters to help with funnel
    prior(student_t(3, 0, 0.5), class = "sd", lb = 0),
    # Less informative correlation prior
    prior(lkj_corr_cholesky(2), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)  # Increased further
)

pairs(fit_hier_age_speed_recommended_v2)

# PICK v3 ----

fit_hier_age_speed_recommended_v3 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # Keep intercept moderately informative
    prior(normal(0, 0.7), class = "Intercept"),
    # Keep informative prior for base velocity
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"),
    # Try regularizing gamma prior for SDs
    prior(gamma(2, 4), class = "sd"),
    # Stronger prior on correlations
    prior(lkj_corr_cholesky(5), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

pairs(fit_hier_age_speed_recommended_v3)

summary(fit_hier_age_speed_recommended_v3)

# PICKED MODELS ----

all_predictors = list()

fit_recommended <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"),
    prior(normal(0.3, 0.3), class = "b", coef = "z_GDS"),
    prior(exponential(3), class = "sd"),
    prior(lkj_corr_cholesky(4), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000
)

fit_hier_age_speed_recommended_v3 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    # Keep intercept moderately informative
    prior(normal(0, 0.7), class = "Intercept"),
    # Keep informative prior for base velocity
    prior(normal(-0.7, 0.3), class = "b", coef = "z_BASE_VELOCITY"),
    # Try regularizing gamma prior for SDs
    prior(gamma(2, 4), class = "sd"),
    # Stronger prior on correlations
    prior(lkj_corr_cholesky(5), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

# Rename existing models
fit_gender_dep <- fit_recommended
fit_age_speed <- fit_hier_age_speed_recommended_v3

# New version of gender+depression model with different priors
fit_gender_dep_v2 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + z_GDS + (1 + z_BASE_VELOCITY + z_GDS | GENDER),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(normal(0, 0.7), class = "Intercept"),
    prior(normal(-0.5, 0.4), class = "b", coef = "z_BASE_VELOCITY"),
    prior(normal(0.2, 0.4), class = "b", coef = "z_GDS"),
    prior(gamma(3, 3), class = "sd"),  # Different shape for variance components
    prior(lkj_corr_cholesky(3), class = "L")  # Less restrictive correlation prior
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)

summary(fit_gender_dep_v2)

# New version of age+speed model with different priors
fit_age_speed_v2 <- brm(
  FALLER ~ 1 + z_BASE_VELOCITY + (1 + z_BASE_VELOCITY | AGE_GROUP),
  data = data,
  family = bernoulli(link = "probit"),
  prior = c(
    prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(-0.6, 0.4), class = "b", coef = "z_BASE_VELOCITY"),
    prior(student_t(3, 0, 0.5), class = "sd", lb = 0),  # Try student_t instead of gamma
    prior(lkj_corr_cholesky(3), class = "L")
  ),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95)
)

summary(fit_age_speed_v2)

all_predictors[["fit_gender_dep"]] <- c("z_BASE_VELOCITY", "z_GDS")
all_predictors[["fit_gender_dep_v2"]] <- c("z_BASE_VELOCITY", "z_GDS")
all_predictors[["fit_age_speed"]] <- c("z_BASE_VELOCITY")
all_predictors[["fit_age_speed_v2"]] <- c("z_BASE_VELOCITY")

# Combine all models into list
fits <- list(
  gender_dep = fit_gender_dep,
  gender_dep_v2 = fit_gender_dep_v2, 
  age_speed = fit_age_speed,
  age_speed_v2 = fit_age_speed_v2
)

get_prior_info <- function(fit, model_name) {
  prior_df <- as.data.frame(prior_summary(fit))
  prior_df$model <- model_name
  return(prior_df)
}

# Apply the function to your list of fits and combine results
prior_df <- do.call(rbind, 
                    Map(get_prior_info, 
                        fits, 
                        names(fits)))

loo_results <- list(
  gender_dep = loo(fit_gender_dep),
  gender_dep_v2 = loo(fit_gender_dep_v2), 
  age_speed = loo(fit_age_speed),
  age_speed_v2 = loo(fit_age_speed_v2)
)

# Extract values
loo_df <- data.frame(
  model = names(loo_results),
  elpd = sapply(loo_results, function(x) x$estimates["elpd_loo", "Estimate"]),
  se = sapply(loo_results, function(x) x$estimates["elpd_loo", "SE"]),
  p_loo = sapply(loo_results, function(x) x$estimates["p_loo", "Estimate"]),
  looic = sapply(loo_results, function(x) x$estimates["looic", "Estimate"])
)

# Extract coefficients and model info
model_stats <- data.frame(
  model = names(fits),
  intercept = sapply(fits, function(x) fixef(x)["Intercept", "Estimate"]),
  intercept_se = sapply(fits, function(x) fixef(x)["Intercept", "Est.Error"]),
  n_obs = sapply(fits, function(x) nobs(x)),
  formula = sapply(fits, function(x) as.character(x$formula)[2])
)

# Extract fixed effects for each model
all_coefs <- do.call(rbind, lapply(names(fits), function(model) {
  coefs <- as.data.frame(fixef(fits[[model]]))
  data.frame(
    model = model,
    parameter = rownames(coefs),
    estimate = coefs$Estimate,
    std_error = coefs$Est.Error
  )
}))

formulas <- list(
  gender_dep = formula(fit_gender_dep),
  gender_dep_v2 = formula(fit_gender_dep_v2), 
  age_speed = formula(fit_age_speed),
  age_speed_v2 = formula(fit_age_speed_v2)
)

# For more readable format:
formulas_df <- data.frame(
  model = names(formulas),
  formula = sapply(formulas, function(x) {
    # Replace long formula with more readable format
    gsub(" \\+ ", "\n+ ", deparse1(x))
  })
)

# Write to CSV
saveRDS(fits, paste0(c(models_path, "faller_classification_models.rds"), collapse = "/"))
write.csv(formulas_df, paste0(c(models_path, "model_coefficients_params.csv"), collapse = "/"), row.names = FALSE)
write.csv(prior_df, paste0(c(models_path, "model_priors.csv"), collapse = "/"), row.names = FALSE)
write.csv(loo_df, paste0(c(results_path, "loo-comparison.csv"), collapse = "/"), row.names = FALSE)
write.csv(all_coefs, paste0(c(results_path, "model_coefficients.csv"), collapse = "/"), row.names = FALSE)
write.csv(model_stats, paste0(c(results_path, "model_summary.csv"), collapse = "/"), row.names = FALSE)

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

# Loop through each fit and create plot
plots_list <- lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  pp_check(fit) + 
    ggtitle(model_name) +
    ylim(0, 2.0) +
    theme_minimal() +
    fits_theme
})

combined_plot <- wrap_plots(plots_list, ncol = 2) + fits_theme
ggsave(paste0(c(plots_path, "all_pp_checks.png"), collapse = "/"), combined_plot, width = 25, height = 20)

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
  mcmc_area_plot <- mcmc_areas(fit) +
    labs(
      title = model_name,
      x = "Standardized Coefficient",
      y = "Parameter",
      subtitle = "Shaded areas represent 95% credible intervals"
    ) +
    base_theme
  mcmc_traces_plot <- mcmc_trace(fit) +
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
  mcmc_intervals_plot <- mcmc_intervals(fit) +
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



# Compare models ----

# Compare models using LOO-CV

comp <- loo_compare(
  loo(fall_risk_models[["gender_dep"]]),
  loo(fall_risk_models[["gender_dep_v2"]]),
  loo(fall_risk_models[["age_speed"]]),
  loo(fall_risk_models[["age_speed_v2"]])
)
print(comp)

# Plot GENDER_DEP comparisons ----

new_data <- expand.grid(
  BASE_VELOCITY = seq(from = min(data$BASE_VELOCITY), 
                      to = max(data$BASE_VELOCITY), 
                      length.out = 100),
  GENDER = c(0, 1),
  GDS = mean(data$GDS)
)

# Create scaled data frame
scaled_new_data <- new_data

# Transform TMT_B if it exists
if ("TMT_B" %in% names(new_data)) {
  scaled_new_data$TMT_B <- log(new_data$TMT_B)
}

# Get columns to scale (excluding GENDER and AGE_GROUP)
cols_to_prefix <- names(new_data)[!names(new_data) %in% c("GENDER", "AGE_GROUP")]

# Scale the selected columns
scaled_new_data[, cols_to_prefix] <- scale(new_data[, cols_to_prefix])

# Add z_ prefix to scaled columns
names(scaled_new_data)[names(scaled_new_data) %in% cols_to_prefix] <- paste0("z_", cols_to_prefix)

compared_fits <- c("gender_dep", "gender_dep_v2")

# Compare posterior predictions
post_preds1 <- posterior_predict(fall_risk_models[[compared_fits[1]]], newdata = scaled_new_data)
post_preds2 <- posterior_predict(fall_risk_models[[compared_fits[2]]], newdata = scaled_new_data)

# Create posterior comparison plot
post_df <- data.frame(
  velocity = rep(new_data$BASE_VELOCITY, 2),
  probability = c(colMeans(post_preds1), colMeans(post_preds2)),
  model = rep(compared_fits, each = nrow(new_data))
)

predict_comp_gender_dep_plot <- ggplot(post_df, aes(x = velocity, y = probability, color = model)) +
  geom_line() +
  ylim(0, 1) +
  labs(x = "Walking Velocity (m/s)",
       y = "Posterior Predicted Probability",
       title = "Posterior Predictions Comparison") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", linewidth = 0)
  )

ggsave(paste0(c(plots_path, paste0("predict_comp_gender_dep.png")), collapse = "/"), predict_comp_gender_dep_plot, width = 10, height = 6)

# Extract and compare key parameters
# Extract parameters correctly
parameters <- data.frame(
  model = rep(compared_fits, each = 4),
  parameter = rep(c("Intercept", "Velocity", "SD_Intercept", "SD_Velocity"), 2),
  estimate = c(
    # Regular model parameters
    fixef(fall_risk_models[[compared_fits[1]]])[1,1],  # Intercept
    fixef(fall_risk_models[[compared_fits[1]]])[2,1],  # Velocity
    as.numeric(VarCorr(fall_risk_models[[compared_fits[1]]])$GENDER$sd[1]),  # SD Intercept
    as.numeric(VarCorr(fall_risk_models[[compared_fits[1]]])$GENDER$sd[2]),  # SD Velocity
    
    # Informative model parameters
    fixef(fall_risk_models[[compared_fits[2]]])[1,1],  # Intercept
    fixef(fall_risk_models[[compared_fits[2]]])[2,1],  # Velocity
    as.numeric(VarCorr(fall_risk_models[[compared_fits[2]]])$GENDER$sd[1]),  # SD Intercept
    as.numeric(VarCorr(fall_risk_models[[compared_fits[2]]])$GENDER$sd[2])   # SD Velocity
  )
)

# Create comparison plot
param_comp_gender_dep_plot <- ggplot(parameters, aes(x = parameter, y = estimate, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Parameter Estimates Comparison",
       x = "Parameter",
       y = "Estimate") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", linewidth = 0)
  )
  coord_flip()

ggsave(paste0(c(plots_path, paste0("param_comp_gender_dep.png")), collapse = "/"), param_comp_gender_dep_plot, width = 10, height = 6)

# Plot AGE_SPEED comparisons ----

new_data <- expand.grid(
  BASE_VELOCITY = seq(from = min(data$BASE_VELOCITY), 
                      to = max(data$BASE_VELOCITY), 
                      length.out = 100),
  AGE_GROUP = 1:3,
  GDS = mean(data$GDS)
)

# Create scaled data frame
scaled_new_data <- new_data

# Transform TMT_B if it exists
if ("TMT_B" %in% names(new_data)) {
  scaled_new_data$TMT_B <- log(new_data$TMT_B)
}

# Get columns to scale (excluding GENDER and AGE_GROUP)
cols_to_prefix <- names(new_data)[!names(new_data) %in% c("GENDER", "AGE_GROUP")]

# Scale the selected columns
scaled_new_data[, cols_to_prefix] <- scale(new_data[, cols_to_prefix])

# Add z_ prefix to scaled columns
names(scaled_new_data)[names(scaled_new_data) %in% cols_to_prefix] <- paste0("z_", cols_to_prefix)

compared_fits <- c("age_speed", "age_speed_v2")

# Compare posterior predictions
post_preds1 <- posterior_predict(fall_risk_models[[compared_fits[1]]], newdata = scaled_new_data)
post_preds2 <- posterior_predict(fall_risk_models[[compared_fits[2]]], newdata = scaled_new_data)

# Create posterior comparison plot
post_df <- data.frame(
  velocity = rep(new_data$BASE_VELOCITY, 2),
  probability = c(colMeans(post_preds1), colMeans(post_preds2)),
  model = rep(compared_fits, each = nrow(new_data))
)

predict_comp_age_speed_plot <- ggplot(post_df, aes(x = velocity, y = probability, color = model)) +
  geom_line() +
  ylim(0, 1) +
  labs(x = "Walking Velocity (m/s)",
       y = "Posterior Predicted Probability",
       title = "Posterior Predictions Comparison") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", linewidth = 0)
  )

ggsave(paste0(c(plots_path, paste0("predict_comp_age_speed.png")), collapse = "/"), predict_comp_age_speed_plot, width = 10, height = 6)

# Extract and compare key parameters
# Extract parameters correctly
parameters <- data.frame(
  model = rep(compared_fits, each = 4),
  parameter = rep(c("Intercept", "Velocity", "SD_Intercept", "SD_Velocity"), 2),
  estimate = c(
    # Regular model parameters
    fixef(fall_risk_models[[compared_fits[1]]])[1,1],  # Intercept
    fixef(fall_risk_models[[compared_fits[1]]])[2,1],  # Velocity
    as.numeric(VarCorr(fall_risk_models[[compared_fits[1]]])$AGE_GROUP$sd[1]),  # SD Intercept
    as.numeric(VarCorr(fall_risk_models[[compared_fits[1]]])$AGE_GROUP$sd[2]),  # SD Velocity
    
    # Informative model parameters
    fixef(fall_risk_models[[compared_fits[2]]])[1,1],  # Intercept
    fixef(fall_risk_models[[compared_fits[2]]])[2,1],  # Velocity
    as.numeric(VarCorr(fall_risk_models[[compared_fits[2]]])$AGE_GROUP$sd[1]),  # SD Intercept
    as.numeric(VarCorr(fall_risk_models[[compared_fits[2]]])$AGE_GROUP$sd[2])   # SD Velocity
  )
)

# Create comparison plot
param_comp_age_speed_plot <- ggplot(parameters, aes(x = parameter, y = estimate, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Parameter Estimates Comparison",
       x = "Parameter",
       y = "Estimate") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", linewidth = 0)
  )
  coord_flip()

ggsave(paste0(c(plots_path, paste0("param_comp_age_speed.png")), collapse = "/"), param_comp_age_speed_plot, width = 10, height = 6)

# Plot AGE_SPEED prediction fit with predicted data points ----

new_data <- expand.grid(
  BASE_VELOCITY = seq(from = min(data$BASE_VELOCITY), 
                      to = max(data$BASE_VELOCITY), 
                      length.out = 100),
  AGE_GROUP = 1:3,
  GDS = mean(data$GDS)
)

# Create scaled data frame
scaled_new_data <- new_data

# Transform TMT_B if it exists
if ("TMT_B" %in% names(new_data)) {
  scaled_new_data$TMT_B <- log(new_data$TMT_B)
}

# Get columns to scale (excluding GENDER and AGE_GROUP)
cols_to_prefix <- names(new_data)[!names(new_data) %in% c("GENDER", "AGE_GROUP")]

# Scale the selected columns
scaled_new_data[, cols_to_prefix] <- scale(new_data[, cols_to_prefix])

# Add z_ prefix to scaled columns
names(scaled_new_data)[names(scaled_new_data) %in% cols_to_prefix] <- paste0("z_", cols_to_prefix)

fit_name <- "age_speed_v2"

fit <- fall_risk_models[[fit_name]]

# Get predictions
preds_epred <- posterior_epred(fit, newdata = scaled_new_data)

predictions <- data.frame(
  new_data,
  estimate = colMeans(preds_epred),
  lower = apply(preds_epred, 2, quantile, probs = 0.025),
  upper = apply(preds_epred, 2, quantile, probs = 0.975)
)

preds_binary <- posterior_predict(fit, newdata = scaled_new_data)
binary_preds <- preds_binary[1,]  # take one draw for binary predictions
predictions$binary = binary_preds

prob_age_speed_plot <- ggplot(predictions, aes(x = BASE_VELOCITY)) +
  # Add probability curves and CI ribbons
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(AGE_GROUP)), alpha = 0.2) +
  geom_line(aes(y = estimate, color = factor(AGE_GROUP))) +
  # Add binary predictions as dots
  geom_dots(data = predictions, 
            aes(y = binary, x = BASE_VELOCITY, color = factor(AGE_GROUP)),
            side = ifelse(predictions$binary == 1, "top", "bottom"),
            dotsize = 1.0,
            pch = 19,
            alpha = 0.5) +
  scale_color_discrete(name = "Age Group",
                       labels = c("< 75", "75-80", "> 80")) +
  scale_fill_discrete(name = "Age Group",
                      labels = c("< 75", "75-80", "> 80")) +
  ylim(0, 1) +
  labs(x = "Walking Velocity (m/s)",
       y = "Probability of being faller", title = paste0(c(fit_name, "Predicted data"), collapse = " | ")) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", linewidth = 0)
  )

ggsave(paste0(c(plots_path, paste0("prop_", fit_name, ".png")), collapse = "/"), prob_age_speed_plot, width = 10, height = 6)


# Plot GENDER_DEP prediction fit with predicted data points ----

new_data <- expand.grid(
  BASE_VELOCITY = seq(from = min(data$BASE_VELOCITY), 
                      to = max(data$BASE_VELOCITY), 
                      length.out = 100),
  GENDER = c(0, 1),
  GDS = mean(data$GDS)
)

# Create scaled data frame
scaled_new_data <- new_data

# Transform TMT_B if it exists
if ("TMT_B" %in% names(new_data)) {
  scaled_new_data$TMT_B <- log(new_data$TMT_B)
}

# Get columns to scale (excluding GENDER and AGE_GROUP)
cols_to_prefix <- names(new_data)[!names(new_data) %in% c("GENDER", "AGE_GROUP")]

# Scale the selected columns
scaled_new_data[, cols_to_prefix] <- scale(new_data[, cols_to_prefix])

# Add z_ prefix to scaled columns
names(scaled_new_data)[names(scaled_new_data) %in% cols_to_prefix] <- paste0("z_", cols_to_prefix)

fit_name <- "gender_dep_v2"

fit <- fall_risk_models[[fit_name]]

# Get predictions
preds_epred <- posterior_epred(fit, newdata = scaled_new_data)

predictions <- data.frame(
  new_data,
  estimate = colMeans(preds_epred),
  lower = apply(preds_epred, 2, quantile, probs = 0.025),
  upper = apply(preds_epred, 2, quantile, probs = 0.975)
)

preds_binary <- posterior_predict(fit, newdata = scaled_new_data)
binary_preds <- preds_binary[1,]  # take one draw for binary predictions
predictions$binary = binary_preds

prob_gender_dep_plot <- ggplot(predictions, aes(x = BASE_VELOCITY)) +
  # Add probability curves and CI ribbons
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(GENDER)), alpha = 0.2) +
  geom_line(aes(y = estimate, color = factor(GENDER))) +
  # Add binary predictions as dots
  geom_dots(data = predictions, 
            aes(y = binary, x = BASE_VELOCITY, color = factor(GENDER)),
            side = ifelse(predictions$binary == 1, "top", "bottom"),
            dotsize = 1.0,
            pch = 19,
            alpha = 0.5) +
  scale_color_discrete(name = "GENDER",
                       labels = c("0 = male", "1 = female")) +
  scale_fill_discrete(name = "GENDER",
                      labels = c("0 = male", "1 = female")) +
  ylim(0, 1) +
  labs(x = "Walking Velocity (m/s)",
       y = "Probability of being faller", title = paste0(c(fit_name, "Predicted data"), collapse = " | ")) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", linewidth = 0)
  )

ggsave(paste0(c(plots_path, paste0("prop_", fit_name, ".png")), collapse = "/"), prob_gender_dep_plot, width = 10, height = 6)
