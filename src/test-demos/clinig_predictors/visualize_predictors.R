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

# Source helper functions ----

source("src/clinical_demog/utils/analysis_plotting_helper_functions.R")
source("src/clinical_demog/utils/analysis_helper_functions.R")


# Load and preview data ----

data <- readRDS("data/clinical_demog_clean.rds")
head(data)
names(data)

# Set common theme ----
fontsize <- 12
common_theme <- theme_minimal() +
  theme(
    text = element_text(size = fontsize),
    axis.title = element_text(size = fontsize),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = fontsize, face = "bold"),
    legend.text = element_text(size = fontsize),
    legend.title = element_text(size = fontsize),
    panel.border = element_rect(color = "grey90", 
                                fill = NA),
    plot.background = element_rect(fill = "white", linewidth = 0),
    panel.spacing.x = unit(0.7, "cm"),
    panel.spacing.y = unit(0.5, "cm"),
    aspect.ratio = 0.7,
    strip.text.y = element_text(angle = 0, 
                                face = "bold", 
                                size = fontsize,
                                margin = margin(r = 10)),
    strip.text.x = element_text(margin = margin(b = 10))
  )


# Catgorize and reshape dataframe for facet plotting ----

# Create the category mapping
# cols_list <- list(
#   DEMO = c("AGE", "GENDER"),
#   BALANCE = c("FSST", "DGI"),
#   GAIT = c("BASE_VELOCITY", "S3_VELOCITY"),
#   COGNITIVE = c("GCS_NEUROTRAX", "TMT_A", "TMT_B"),
#   DEPRESSION = c("GDS")
# )
cols_list <- list(
  PHYSICAL = c("GENDER", "DGI", "FSST", "BASE_VELOCITY", "S3_VELOCITY"),
  COGNITIVE = c("GCS_NEUROTRAX", "TMT_A", "TMT_B"),
  DEPRESSION = c("GDS")
)
predictor_category_names <- names(cols_list)

predictor_categories <- rep(names(cols_list), times = sapply(cols_list, length))
names(predictor_categories) <- unlist(cols_list)
predictors <- names(predictor_categories)


df_long <- data %>%
  pivot_longer(
    cols = all_of(names(data)[!names(data) %in% c("ID", "FALLER", "DATE_OF_EVALUATION")]),
    names_to = "PREDICTOR",
    values_to = "PREDICTOR_VALUE"
  ) %>%
  select(ID, PREDICTOR, PREDICTOR_VALUE, FALLER) %>%
  mutate(PREDICTOR_CATEGORY = predictor_categories[PREDICTOR])

plot_data <- df_long[!is.na(df_long$PREDICTOR_CATEGORY), ] %>%
  mutate(
    PREDICTOR_CATEGORY = factor(PREDICTOR_CATEGORY, 
                                levels = predictor_category_names),
    PREDICTOR = factor(PREDICTOR,
                       levels = predictors),
    FALLER = factor(FALLER)
  )

ggplot(data = plot_data) +
  geom_point(mapping = aes(x = PREDICTOR_VALUE, y = FALLER, colour = FALLER), size = 3, show.legend = TRUE, alpha = 1/10) +
  labs(title='Predictors by Category') +
  facet_grid(PREDICTOR_CATEGORY ~ PREDICTOR, 
             switch = 'y', scales = 'free_x') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  common_theme
  
ggsave('plots/predictors/predictors.png', width = 20, height = 5)
# Basic distribution checks for TMT_A and TMT_B ----

source("src/test-demos/clinig_predictors/visualize_predictors_helper_functions.R")

# Run the function

compare_log_transformed(data, "TMT_A")
ggsave('plots/predictors/TMT_A.png', width = 20, height = 20)
compare_log_transformed(data, "TMT_B")
ggsave('plots/predictors/TMT_B.png', width = 20, height = 20)


# Create log transformed variables  ----

# Create the category mapping
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
# First scale the data as you did
data_scaled <- data
data_scaled[paste0("z_", predictors[predictors != "GENDER"])] <- scale(data[predictors[predictors != "GENDER"]])

# Create category mapping to take into account scaled predictors
cols_list <- list(
  PHYSICAL = c("GENDER", paste0("z_", c("DGI", "FSST", "BASE_VELOCITY", "S3_VELOCITY"))),
  COGNITIVE = paste0("z_", c("GCS_NEUROTRAX", "log_TMT_A", "log_TMT_B")),
  DEPRESSION = paste0("z_", c("GDS"))
)
predictor_category_names <- names(cols_list)

predictor_categories <- rep(names(cols_list), times = sapply(cols_list, length))
names(predictor_categories) <- unlist(cols_list)
predictors <- names(predictor_categories)


df_long <- data_scaled %>%
  pivot_longer(
    cols = all_of(names(data_scaled)[!names(data_scaled) %in% c("ID", "FALLER", "DATE_OF_EVALUATION")]),
    names_to = "PREDICTOR",
    values_to = "PREDICTOR_VALUE"
  ) %>%
  select(ID, PREDICTOR, PREDICTOR_VALUE, FALLER) %>%
  mutate(PREDICTOR_CATEGORY = predictor_categories[PREDICTOR])

plot_data <- df_long[!is.na(df_long$PREDICTOR_CATEGORY), ] %>%
  mutate(
    PREDICTOR_CATEGORY = factor(PREDICTOR_CATEGORY, 
                                levels = predictor_category_names),
    PREDICTOR = factor(PREDICTOR,
                       levels = predictors),
    FALLER = factor(FALLER)
  )

ggplot(data = plot_data) +
  geom_point(mapping = aes(x = PREDICTOR_VALUE, y = FALLER, colour = FALLER), size = 3, show.legend = TRUE, alpha = 1/10) +
  labs(title='Z-scored Predictors by Category') +
  facet_grid(PREDICTOR_CATEGORY ~ PREDICTOR, 
             switch = 'y') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  common_theme
ggsave('plots/predictors/z-scored_predictors.png', width = 20, height = 5)


# Bayesian areas ----

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

data_scaled <- data
# data_scaled[paste0("z_", predictors[predictors != "GENDER"])] <- scale(data[predictors[predictors != "GENDER"]])
data_scaled[paste0("z_", predictors)] <- scale(data[predictors])

cols_list <- list(
  PHYSICAL = paste0("z_", c("GENDER", "DGI", "FSST", "BASE_VELOCITY", "S3_VELOCITY")),
  COGNITIVE = paste0("z_", c("GCS_NEUROTRAX", "log_TMT_A", "log_TMT_B")),
  DEPRESSION = paste0("z_", c("GDS"))
)
predictor_category_names <- names(cols_list)

predictor_categories <- rep(names(cols_list), times = sapply(cols_list, length))
names(predictor_categories) <- unlist(cols_list)
predictors <- names(predictor_categories)

z_data <- data_scaled

z_data$FALLER <- factor(z_data$FALLER)

fall_class_formula <- bf(as.formula(paste("FALLER ~ 1 +", paste(predictors, collapse = " + "))), family = "bernoulli")

# Example with specific priors
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
  data = z_data,
  prior = fall_class_priors
)

summary(fall_class_fit)


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

loo_result <- loo(fall_class_fit)
print(loo_result)