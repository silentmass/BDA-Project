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

default_prior_b <- prior(student_t(3, 0, 2.5), class = "b")

fit_prefix <- "faller_classification_student-t-prior"

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

# Scaled predictor names
cols_list <- list(
  PHYSICAL = paste0("z_", original_cols_list[["PHYSICAL"]]),
  SPEED = paste0("z_", original_cols_list[["SPEED"]]),
  COGNITIVE = paste0("z_", original_cols_list[["COGNITIVE"]]),
  DEPRESSION = paste0("z_", original_cols_list[["DEPRESSION"]])
)

predictor_categories <- get_predictor_categories_from_cols_list(cols_list)

# Init lists for formulas, priors, and fits
all_predictors <- list()
formulas <- list()
priors <- list()
fits <- list()
summaries <- list()
loo_results <- list()


# Training models ----
# Running this section trains all the models
# Remember to save results separately

#### All categories ----

selected_categories <- names(cols_list)
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(all_predictors[[fit_name]], collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Gait speed - stronger prior based on your data showing clear difference
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  default_prior_b
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

#### Physical and Speed ----

selected_categories <- c("PHYSICAL", "SPEED")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(all_predictors[[fit_name]], collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Gait speed - stronger prior based on your data showing clear difference
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  default_prior_b
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


#----

#### Physical ----

selected_categories <- c("PHYSICAL")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(all_predictors[[fit_name]], collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  default_prior_b
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


#### Depression----

selected_categories <- c("DEPRESSION")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(all_predictors[[fit_name]], collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  default_prior_b
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


#### Speed ----

selected_categories <- c("SPEED")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(all_predictors[[fit_name]], collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Gait speed - stronger prior based on your data showing clear difference
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept")
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


#### Cognitive ----

selected_categories <- c("COGNITIVE")
fit_name <- paste0("c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(paste("FALLER ~ 1 +", paste(all_predictors[[fit_name]], collapse = " + "))), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  default_prior_b
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


#### Cognitive with spline ----

selected_categories <- c("COGNITIVE")
fit_name <- paste0("spline-c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(
  paste("FALLER ~ 1 +", 
        paste(all_predictors[[fit_name]], collapse = " + "), 
        paste0(" + s(", paste(all_predictors[[fit_name]], collapse = ") + s("), ")"))
), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  default_prior_b
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


#### Hierarchical by age ----

hierarchical_col <- "AGE_GROUP"
selected_categories <- names(cols_list)
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(
  paste("FALLER ~ 1 +", 
        paste(all_predictors[[fit_name]], collapse = " + "), 
        "+ (1 +", paste(all_predictors[[fit_name]], collapse = " + "), paste0("| ", hierarchical_col, ")"))
), family = "bernoulli")

priors[[fit_name]] <- c(
  # Specific priors first
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "sd", group = "AGE_GROUP", coef = "z_BASE_VELOCITY"),
  # Defaults last
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd"),
  default_prior_b
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


#### Hierarchical by gender ----

hierarchical_col <- "GENDER"
selected_categories <- names(cols_list)
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(
  paste("FALLER ~ 1 +", 
        paste(all_predictors[[fit_name]], collapse = " + "), 
        "+ (1 +", paste(all_predictors[[fit_name]], collapse = " + "), paste0("| ", hierarchical_col, ")"))
), family = "bernoulli")

priors[[fit_name]] <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "sd", group = "GENDER", coef = "z_BASE_VELOCITY"),
  # Other variables - more uncertain
  # Defaults last
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd"),
  default_prior_b
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


#### Hierarchical by gender, only speed and depression ----

hierarchical_col <- "GENDER"
selected_categories <- c("SPEED", "DEPRESSION")
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(
  paste("FALLER ~ 1 +",
        paste(all_predictors[[fit_name]], collapse = " + "),
        "+ (1 +", paste(all_predictors[[fit_name]], collapse = " + "), paste0("| ", hierarchical_col, ")"))
), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  # Defaults last
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1.5), class = "sd"),
  default_prior_b
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


#### Hierarchical model using only speed as parameter, age group hierarchy ----

hierarchical_col <- "AGE_GROUP"
selected_categories <- c("SPEED")
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(
  paste("FALLER ~ 1 +", 
        paste(all_predictors[[fit_name]], collapse = " + "), 
        "+ (1 +", paste(all_predictors[[fit_name]], collapse = " + "), paste0("| ", hierarchical_col, ")"))
), family = "bernoulli")

priors[[fit_name]] <- c(
  prior(normal(-0.5, 0.5), class = "sd", group = "AGE_GROUP", coef = "z_BASE_VELOCITY"),
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


#### Hierarchical by age, only cognition category, using spline ----

hierarchical_col <- "AGE_GROUP"
selected_categories <- c("COGNITIVE")
fit_name <- paste0(paste0("hierarchical-", hierarchical_col), "-spline-c-", paste0(selected_categories, collapse = "_"))

all_predictors[[fit_name]] <- unique(names(predictor_categories[predictor_categories %in% selected_categories]))

formulas[[fit_name]] <- bf(as.formula(
  paste("FALLER ~ 1 +",
        paste0("s(", paste(all_predictors[[fit_name]], collapse = ") + s("), ")"),
        "+ (1 + ", paste(all_predictors[[fit_name]], collapse = " + "), paste0("| ", hierarchical_col, ")"))
), family = "bernoulli")

priors[[fit_name]] <- c(
  # Other variables - more uncertain
  prior(normal(0, 1), class = "Intercept"),
  default_prior_b
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

# For more readable format:
formulas_df <- data.frame(
  model = names(formulas),
  formula = sapply(formulas, function(x) {
    # Replace long formula with more readable format
    gsub(" \\+ ", "\n+ ", deparse1(x))
  })
)

# Write to CSV
write.csv(loo_df, paste0(c(results_path, "loo-comparison.csv"), collapse = "/"), row.names = FALSE)
write.csv(all_coefs, paste0(c(results_path, "model_coefficients.csv"), collapse = "/"), row.names = FALSE)
write.csv(model_stats, paste0(c(results_path, "model_summary.csv"), collapse = "/"), row.names = FALSE)
saveRDS(fits, paste0(c(models_path, "faller_classification_models.rds"), collapse = "/"))
write.csv(formulas_df, paste0(c(models_path, "model_coefficients_params.csv"), collapse = "/"), row.names = FALSE)


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
  pp_check(fit) + 
    ggtitle(format_model_name(model_name)) +
    ylim(0, 2.0) +
    theme_minimal() +
    fits_theme
})

combined_plot <- wrap_plots(plots_list, ncol = 2) + fits_theme
ggsave(paste0(c(plots_path, "all_pp_checks.png"), collapse = "/"), combined_plot, width = 25, height = 20)


# Loop through each fit and create plot
plots_list <- lapply(names(fits), function(model_name) {
  fit <- fits[[model_name]]
  mcmc_area_plot <- plot_mcmc(fit, model_name, all_predictors[[model_name]], "area")
  mcmc_traces_plot <- plot_mcmc(fit, model_name, all_predictors[[model_name]], "trace")
  mcmc_intervals_plot <- plot_mcmc(fit, model_name, all_predictors[[model_name]], "interval")
  
  ggsave(paste0(c(mcmc_path, "areas", paste0("areas_", model_name, ".png")), collapse = "/"), mcmc_area_plot, width = 20, height = 20)
  ggsave(paste0(c(mcmc_path, "traces", paste0("traces_", model_name, ".png")), collapse = "/"), mcmc_traces_plot, width = 20, height = 20)
  ggsave(paste0(c(mcmc_path, "intervals", paste0("intervals_", model_name, ".png")), collapse = "/"), mcmc_intervals_plot, width = 20, height = 20)
})

# projpred variable selection ----

#here validate_search=TRUE, even if it is slower that way. use fit from all parameters, normal(0,1) prior
#because the results looked better (elpd_loo was best with that prior)
fitmodel = fits[["c-PHYSICAL_SPEED_COGNITIVE_DEPRESSION"]]
summary(fitmodel)
varsel <- cv_varsel(fitmodel, method='forward', cv_method='loo', validate_search=TRUE)
plot(varsel, stats = c('elpd', 'pctcorr'), deltas=FALSE, text_angle = 45)
ggsave("plots/fall_variable_selection_prior_0_1.png")
# ggsave("plots/fall_variable_selection_prior_0_5.png")
nsel<-suggest_size(varsel)
#-------- 

