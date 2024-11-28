# Load packages and set theme ----

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

# Set SEED ----
SEED = 2024

# Source helper functions ----

source("src/clinical_demog/utils/analysis_plotting_helper_functions.R")
source("src/clinical_demog/utils/analysis_helper_functions.R")

# Load and preview data ----

data <- readRDS("data/clinical_demog_clean.rds")
head(data)
names(data)

data$FALLER <- factor(data$FALLER)

data <- data %>%
  mutate(AGE_GROUP = case_when(
    AGE < 75 ~ 1,
    AGE >= 75 & AGE <= 80 ~ 2,
    AGE > 80 ~ 3
  ))

# Log transform TMT and z-score scale and fit ----

physical_cols <- c("AGE", "GENDER", "DGI", "TUG", "FSST")
speed_cols <- c("BASE_VELOCITY", "S3_VELOCITY")
cognitive_cols <- c("GCS_NEUROTRAX", "log_TMT_A", "log_TMT_B")
depression_cols <- c("GDS")

original_cols_list <- list(
  PHYSICAL = physical_cols,
  SPEED = speed_cols,
  COGNITIVE = cognitive_cols,
  DEPRESSION = depression_cols
)

original_predictor_categories <- get_predictor_categories_from_cols_list(original_cols_list)
original_predictors <- names(original_predictor_categories)

# Log transform TMTs
data$log_TMT_A <- log(data$TMT_A)
data$log_TMT_B <- log(data$TMT_B)

# Scale all predictors including log-transformed TMT variables
data[paste0("z_", original_predictors)] <- scale(data[original_predictors]) # Scale also GENDER

# Scaled predictor names
cols_list <- list(
  PHYSICAL = paste0("z_", physical_cols),
  SPEED = paste0("z_", speed_cols),
  COGNITIVE = paste0("z_", cognitive_cols),
  DEPRESSION = paste0("z_", depression_cols)
)

predictor_categories <- get_predictor_categories_from_cols_list(cols_list)

# Init lists for formulas, priors, and fits
all_predictors <- list()
formulas <- list()
priors <- list()
fits <- list()
summaries <- list()
loo_results <- list()

# All categories ----

selected_categories <- names(cols_list)
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Gait speed - stronger prior based on your data showing clear difference
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "z_S3_VELOCITY"),
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# z_FSST only because it was the one selected by cv_varsel ----

fit_name <- "FSST"
predictors <- c("z_FSST")

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])

# Physical ----

selected_categories <- c("PHYSICAL")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])

# Depression----

selected_categories <- c("DEPRESSION")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Speed ----

selected_categories <- c("SPEED")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Gait speed - stronger prior based on your data showing clear difference
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "z_S3_VELOCITY"),
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Cognitive ----

selected_categories <- c("COGNITIVE")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Cognitive with spline ----

selected_categories <- c("COGNITIVE")
fit_name <- paste0("spline-c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "), paste0(" + s(", paste(predictors, collapse = ") + s("), ")"))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Hierarchical by age ----

hierarchical_col <- "AGE_GROUP"
selected_categories <- names(cols_list)
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "), "+ (1 +", paste(predictors, collapse = " + "), paste0("| ", hierarchical_col, ")"))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Specific priors first
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "z_S3_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "sd", group = "AGE_GROUP", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "sd", group = "AGE_GROUP", coef = "z_S3_VELOCITY"),
  # Defaults last
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Hierarchical by gender ----

hierarchical_col <- "GENDER"
selected_categories <- names(cols_list)
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "), "+ (1 +", paste(predictors, collapse = " + "), paste0("| ", hierarchical_col, ")"))), family = "bernoulli")

priors[[fit_name]] <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "z_S3_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "sd", group = "GENDER", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "sd", group = "GENDER", coef = "z_S3_VELOCITY"),
  # Other variables - more uncertain
  # Defaults last
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Hierarchical by gender, only speed and depression ----

hierarchical_col <- "GENDER"
selected_categories <- c("SPEED", "DEPRESSION")
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "), "+ (1 +", paste(predictors, collapse = " + "), paste0("| ", hierarchical_col, ")"))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  # Defaults last
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1.5), class = "sd"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])
save_fit_summary(summaries[[fit_name]], fit_name = fit_name)

loo_results[[fit_name]] <- loo(fits[[fit_name]])
save_loo_result(loo_results[[fit_name]], fit_name = fit_name)

print(fit_name)
print(loo_results[[fit_name]])


# Hierarchical model using only speed as parameter, age group hierarchy ----

hierarchical_col <- "AGE_GROUP"
selected_categories <- c("SPEED")
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "), "+ (1 +", paste(predictors, collapse = " + "), paste0("| ", hierarchical_col, ")"))), family = "bernoulli")

priors[[fit_name]] <- c(
  prior(normal(-0.5, 0.5), class = "sd", group = "AGE_GROUP", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "sd", group = "AGE_GROUP", coef = "z_S3_VELOCITY"),
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1.5), class = "sd")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Hierarchical by age, only cognition category, using spline ----

hierarchical_col <- "AGE_GROUP"
selected_categories <- c("COGNITIVE")
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-spline-c-", paste0(selected_categories, collapse = "_"))

predictors <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

all_predictors[[fit_name]] <- predictors

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste0("s(", paste(predictors, collapse = ") + s("), ")"), "+ (1 + ", paste(predictors, collapse = " + "), paste0("| ", hierarchical_col, ")"))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b")
)

fits[[fit_name]] <- brm(
  formula = formulas[[fit_name]],
  data = data,
  prior = priors[[fit_name]],
  seed = SEED
)

summaries[[fit_name]] <- summary(fits[[fit_name]])

loo_results[[fit_name]] <- loo(fits[[fit_name]])

print(fit_name)
print(loo_results[[fit_name]])


# Combine result into df and save ----

# Extract values
loo_df <- data.frame(
  model = names(loo_results),
  elpd = sapply(loo_results, function(x) x$estimates["elpd_loo", "Estimate"]),
  se = sapply(loo_results, function(x) x$estimates["elpd_loo", "SE"]),
  p_loo = sapply(loo_results, function(x) x$estimates["p_loo", "Estimate"]),
  looic = sapply(loo_results, function(x) x$estimates["looic", "Estimate"])
)

# Export
write.csv(loo_df, paste0(c("results", "loo-comparison.csv"), collapse = "/"), row.names = FALSE)

# Save fits

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

# Write to CSV
write.csv(all_coefs, paste0(c("results", "model_coefficients.csv"), collapse = "/"), row.names = FALSE)
write.csv(model_stats, paste0(c("results", "model_summary.csv"), collapse = "/"), row.names = FALSE)

# Plot fits


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
  # plot_mcmc(fit, model_name, all_predictors[[model_name]])
  pp_check(fit) + 
    ggtitle(gsub("hierarchical", "h", model_name)) +
    ylim(0, 2.0) +
    theme_minimal() +
    fits_theme
})

combined_plot <- wrap_plots(plots_list, ncol = 2) + fits_theme
ggsave(paste0(c("results", "all_pp_checks.png"), collapse = "/"), combined_plot, width = 20, height = 20)


# projpred variable selection ----

#here validate_search=TRUE, even if it is slower that way. use fall_class_fit (not horseshoe) 
#because the results looked better (elpd_loo was smaller (-50) than with horseshoe (-49))
# varsel2 <- cv_varsel(fall_class_fit, method='forward', cv_method='loo', validate_search=FALSE)
# plot(varsel2, stats = c('elpd', 'pctcorr'), deltas=FALSE, text_angle = 45)
# ggsave("plots/fall_variable_selection_prior_0_1.png")
# ggsave("plots/fall_variable_selection_prior_0_5.png")
# nsel<-suggest_size(varsel2)
# (vsel<-solution_terms(varsel2)[1:nsel])
# proj2 <- project(varsel2, nv = nsel, ns = 4000)
# proj2draws <- as.matrix(proj2)
# colnames(proj2draws) <- c("Intercept",vsel)
# round(colMeans(proj2draws),1)
# round(posterior_interval(proj2draws),1)
# 
# mcmc_areas(proj2draws, prob = 0.95, prob_outer = 1,
#            pars = c('Intercept', vsel))

# Plot pp check for all ----

# plot_pp_check(fall_class_fit, names(predictor_categories))
# plot_pp_check(fall_class_horseshoe_fit, names(predictor_categories))
# plot_pp_check(fall_class_fsst_fit, "z_FSST")
# plot_pp_check(fall_class_physical_fit, cols_list$PHYSICAL)
# plot_pp_check(fall_class_speed_fit, cols_list$SPEED)
# plot_pp_check(fall_class_cognition_fit, cols_list$COGNITION)
# plot_pp_check(fall_class_depression_fit, cols_list$DEPRESSION)
# plot_pp_check(fall_class_hierarchical_fit, predictors_for_hier)
# plot_pp_check(fall_class_hierarchical_depr_fit, predictors_depr)

# Save results ----

# Extract fixed effects
# fixed_effects <- fixef(fall_class_fit)
# fixed_effects_df <- as.data.frame(fixed_effects)
# 
# # Add row names as a column
# fixed_effects_df$Parameter <- rownames(fixed_effects_df)
# rownames(fixed_effects_df) <- NULL
# 
# # Save to CSV
# write.csv(fixed_effects_df, "results/fall_model_fixed_effects.csv", row.names = FALSE)

# If you want to save random effects (if any):
# random_effects <- ranef(fall_class_fit)
# write.csv(as.data.frame(random_effects), "fall_model_random_effects.csv")

# To save the full posterior samples:
# posterior_samples <- as.data.frame(posterior)
# write.csv(posterior_samples, "results/fall_model_posterior_samples.csv", row.names = FALSE)

# To save model fit statistics
# fit_stats <- data.frame(
#   elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
#   se_elpd_loo = loo_result$estimates["elpd_loo", "SE"],
#   p_loo = loo_result$estimates["p_loo", "Estimate"],
#   looic = loo_result$estimates["looic", "Estimate"]
# )
# write.csv(fit_stats, "results/fall_model_fit_statistics.csv", row.names = FALSE)
