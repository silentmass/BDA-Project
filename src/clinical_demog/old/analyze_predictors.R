######################### Load packages and set theme ----

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

# if (!require(ggplot2)) {
#   install.packages("ggplot2")
#   library(ggplot2)
# }

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

cmdstan_installed <- function(){
  res <- try(out <- cmdstanr::cmdstan_path(), silent = TRUE)
  !inherits(res, "try-error")
}

if(!cmdstan_installed()){
  install_cmdstan()
}

ggplot2::theme_set(ggplot2::theme_minimal())

#########################  Source helper functions ----

source("src/clinical_demog/utils/analysis_plotting_helper_functions.R")
source("src/clinical_demog/utils/analysis_helper_functions.R")


#########################  Load and preview data ----

data <- readRDS("data/clinical_demog_clean.rds")
head(data)
names(data)

######################### Set common theme ----

mcmc_fontsize <- 12
mcmc_theme <- theme_minimal() +
  theme(
    text = element_text(size = mcmc_fontsize),
    axis.title = element_text(size = mcmc_fontsize),
    axis.text = element_text(size = mcmc_fontsize),
    strip.text = element_text(size = mcmc_fontsize, face = "bold"),
    legend.text = element_text(size = mcmc_fontsize),
    legend.title = element_text(size = mcmc_fontsize),
    panel.border = element_rect(color = "grey90", 
                                fill = NA),
    plot.background = element_rect(fill = "white", linewidth = 0),
    panel.spacing.x = unit(0.7, "cm"),
    panel.spacing.y = unit(0.5, "cm"),
    aspect.ratio = 0.7,
    strip.text.y = element_text(angle = 0, 
                                face = "bold", 
                                size = mcmc_fontsize,
                                margin = margin(r = 10)),
    strip.text.x = element_text(margin = margin(b = 10))
  )

#########################  Log transform TMT and z-score scale and fit ----

cols_list <- list(
  PHYSICAL = c("GENDER", "DGI", "FSST", "BASE_VELOCITY", "S3_VELOCITY"),
  COGNITIVE = c("GCS_NEUROTRAX", "log_TMT_A", "log_TMT_B"),
  DEPRESSION = c("GDS")
)
predictor_category_names <- names(cols_list)

predictor_categories <- rep(names(cols_list), times = sapply(cols_list, length))
names(predictor_categories) <- unlist(cols_list)
predictors <- names(predictor_categories)

# Log transform TMTs
data$log_TMT_A <- log(data$TMT_A)
data$log_TMT_B <- log(data$TMT_B)

# Scale all predictors including log-transformed TMT variables
data[paste0("z_", predictors)] <- scale(data[predictors]) # Scale also GENDER

# Scaled predictor names
cols_list <- list(
  PHYSICAL = paste0("z_", c("GENDER", "DGI", "FSST", "BASE_VELOCITY", "S3_VELOCITY")),
  COGNITIVE = paste0("z_", c("GCS_NEUROTRAX", "log_TMT_A", "log_TMT_B")),
  DEPRESSION = paste0("z_", c("GDS"))
)
predictor_category_names <- names(cols_list)

predictor_categories <- rep(names(cols_list), times = sapply(cols_list, length))
names(predictor_categories) <- unlist(cols_list)
predictors <- names(predictor_categories)

data$FALLER <- factor(data$FALLER)

fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

fall_class_priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  # Gait speed - stronger prior based on your data showing clear difference
  prior(normal(-0.5, 0.5), class = "b", coef = "z_BASE_VELOCITY"),
  prior(normal(-0.5, 0.5), class = "b", coef = "z_S3_VELOCITY"),
  # Other variables - more uncertain
  prior(normal(0, 1), class = "b")
)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = data,
  prior = fall_class_priors
)

summary(fall_class_fit)

posterior <- as_draws_df(fall_class_fit)

mcmc_areas(
  posterior,
  pars = paste0("b_", predictors)
) + mcmc_theme
ggsave("plots/predictors/mcmc_areas.png")

mcmc_trace(
  posterior,
  pars = paste0("b_", predictors)
) + mcmc_theme
ggsave("plots/predictors/mcmc_traces.png")

mcmc_intervals(
  posterior,
  pars = paste0("b_", predictors),
) + mcmc_theme
ggsave("plots/predictors/mcmc_intervals.png")

#projpred variable selection----
varsel2 <- cv_varsel(fall_class_fit, method='forward', cv_method='loo', validate_search=FALSE)
plot(varsel2, stats = c('elpd', 'pctcorr'), deltas=FALSE, text_angle = 45)
#----
loo_result <- loo(fall_class_fit)
print(loo_result)

# Save results

# Extract fixed effects
fixed_effects <- fixef(fall_class_fit)
fixed_effects_df <- as.data.frame(fixed_effects)

# Add row names as a column
fixed_effects_df$Parameter <- rownames(fixed_effects_df)
rownames(fixed_effects_df) <- NULL

# Save to CSV
write.csv(fixed_effects_df, "results/fall_model_fixed_effects.csv", row.names = FALSE)

# If you want to save random effects (if any):
# random_effects <- ranef(fall_class_fit)
# write.csv(as.data.frame(random_effects), "fall_model_random_effects.csv")

# To save the full posterior samples:
posterior_samples <- as.data.frame(posterior)
write.csv(posterior_samples, "results/fall_model_posterior_samples.csv", row.names = FALSE)

# To save model fit statistics
fit_stats <- data.frame(
  elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
  se_elpd_loo = loo_result$estimates["elpd_loo", "SE"],
  p_loo = loo_result$estimates["p_loo", "Estimate"],
  looic = loo_result$estimates["looic", "Estimate"]
)
write.csv(fit_stats, "results/fall_model_fit_statistics.csv", row.names = FALSE)
