####### Run full ----

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
source("src/clinical_demog/utils/visualize_predictors_helper_functions.R")

# Load and preview data ----

data <- readRDS("data/clinical_demog_clean.rds")
head(data)
names(data)

# Set common theme and plot -----

width <- 20
height <- 6

fontsize <- 11
common_theme <- theme_minimal() +
  theme(
    text = element_text(size = fontsize),
    axis.title = element_text(size = fontsize),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = fontsize, face = "bold"),
    legend.text = element_text(size = fontsize),
    legend.title = element_text(size = fontsize),
    panel.border = element_rect(color = "grey90", 
                                fill = NA),
    plot.background = element_rect(fill = "white", linewidth = 0),
    panel.spacing.x = unit(0.6, "cm"),
    panel.spacing.y = unit(0.3, "cm"),
    aspect.ratio = 0.6,
    strip.text.y = element_text(angle = 0, 
                                face = "bold", 
                                size = fontsize,
                                margin = margin(r = 10)),
    strip.text.x = element_text(margin = margin(b = 10))
  )


# Categorize and reshape dataframe for facet plotting

# Set original predictor categories
physical_cols <- c("AGE", "GENDER", "DGI", "TUG", "FSST")
speed_cols <- c("BASE_VELOCITY")
cognitive_cols <- c("GCS_NEUROTRAX", "TMT_B")
depression_cols <- c("GDS") 

original_cols_list <- list(
  PHYSICAL = physical_cols,
  SPEED = speed_cols,
  COGNITIVE = cognitive_cols,
  DEPRESSION = depression_cols
)

original_predictor_categories <- get_predictor_categories_from_cols_list(original_cols_list)
original_predictors <- names(original_predictor_categories)


df_long <- data %>%
  pivot_longer(
    cols = all_of(names(data)[!names(data) %in% c("ID", "FALLER", "DATE_OF_EVALUATION")]),
    names_to = "PREDICTOR",
    values_to = "PREDICTOR_VALUE"
  ) %>%
  select(ID, PREDICTOR, PREDICTOR_VALUE, FALLER) %>%
  mutate(PREDICTOR_CATEGORY = original_predictor_categories[PREDICTOR])

plot_data <- df_long[!is.na(df_long$PREDICTOR_CATEGORY), ] %>%
  mutate(
    PREDICTOR_CATEGORY = factor(PREDICTOR_CATEGORY, 
                                levels = names(original_cols_list)),
    PREDICTOR = factor(PREDICTOR,
                       levels = original_predictors),
    FALLER = factor(FALLER)
  )

# Calculate medians per category and predictor
medians_df <- plot_data %>%
  group_by(PREDICTOR_CATEGORY, PREDICTOR, FALLER) %>%
  summarise(median_val = median(PREDICTOR_VALUE, na.rm = TRUE))

# Add to plot
ggplot(data = plot_data) +
  geom_point(mapping = aes(x = PREDICTOR_VALUE, y = FALLER, colour = FALLER), 
             size = 3, alpha = 1/10) +
  geom_vline(data = medians_df,
             aes(xintercept = median_val, color = FALLER),
             linetype = "dashed") +
  labs(title='Predictors by Category') +
  facet_grid(PREDICTOR_CATEGORY ~ PREDICTOR, 
             switch = 'y', scales = 'free_x') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  common_theme
  
ggsave('plots/predictors/predictors_raw.png', width = width, height = height)

# Basic distribution checks for TMT_A and TMT_B
compare_log_transformed(data, "TMT_A")
ggsave('plots/predictors/TMT_A.png', width = 20, height = 20)
compare_log_transformed(data, "TMT_B")
ggsave('plots/predictors/TMT_B.png', width = 20, height = 20)


# Log transform and z-score scale and fit

transformed_predictor_categories <- get_predictor_categories_from_cols_list(list(
  PHYSICAL = physical_cols,
  SPEED = speed_cols,
  COGNITIVE = sub("^(TMT_[AB])", "log_\\1", cognitive_cols),
  DEPRESSION = depression_cols
))

# Log transform TMTs
data$log_TMT_A <- log(data$TMT_A)
data$log_TMT_B <- log(data$TMT_B)

# Scale all predictors including log-transformed TMT variables
# First scale the data as you did
data_scaled <- data
data_scaled[paste0("z_", names(transformed_predictor_categories))] <- scale(data[names(transformed_predictor_categories)])

# Create category mapping to take into account scaled predictors
cols_list <- list(
  PHYSICAL = paste0("z_", physical_cols),
  SPEED = paste0("z_", speed_cols),
  COGNITIVE = paste0("z_", sub("^(TMT_[AB])", "log_\\1", cognitive_cols)),
  DEPRESSION = paste0("z_", depression_cols)
)


predictor_categories <- get_predictor_categories_from_cols_list(cols_list)


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
    PREDICTOR_CATEGORY = factor(PREDICTOR_CATEGORY, levels = names(cols_list)),
    PREDICTOR = factor(PREDICTOR, levels = names(predictor_categories)),
    FALLER = factor(FALLER)
  )

medians_df <- plot_data %>%
  group_by(PREDICTOR_CATEGORY, PREDICTOR, FALLER) %>%
  summarise(median_val = median(PREDICTOR_VALUE, na.rm = TRUE))

ggplot(data = plot_data) +
  geom_point(mapping = aes(x = PREDICTOR_VALUE, y = FALLER, colour = FALLER), size = 3, show.legend = TRUE, alpha = 1/10) +
  labs(title='Z-scored Predictors by Category') +
  geom_vline(data = medians_df,
             aes(xintercept = median_val, color = FALLER),
             linetype = "dashed") +
  facet_grid(PREDICTOR_CATEGORY ~ PREDICTOR, 
             switch = 'y') +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  common_theme
ggsave('plots/predictors/z-scored_predictors.png', width = width, height = height)
