
####### Run pre ----

#### Load models and libraries ----

# prediction_script.R
library(brms)
library(tidyverse)
library(pROC)

# Set SEED ----
SEED = 2024

# Source helper functions ----
utils_path <- paste0(c("src", "clinical_demog", "utils"), collapse = "/")

source(paste0(c(utils_path, "analysis_plotting_helper_functions.R"), collapse = "/"))
source(paste0(c(utils_path, "analysis_helper_functions.R"), collapse = "/"))

# Load models
fits <- readRDS("models/faller_classification_normal-prior_0-1/faller_classification_models.rds")

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

fit <- fits[["c-SPEED"]]

# Create prediction grid
new_data <- expand.grid(
  z_BASE_VELOCITY = seq(from = min(data$z_BASE_VELOCITY), 
                    to = max(data$z_BASE_VELOCITY), 
                    length.out = 100)
)

# Get posterior predictions of probabilities
preds <- posterior_epred(fit, newdata = new_data)

# Calculate summary statistics
predictions <- data.frame(
  new_data,
  # Mean prediction (probability)
  estimate = colMeans(preds),
  # Uncertainty intervals
  lower = apply(preds, 2, quantile, probs = 0.025),
  upper = apply(preds, 2, quantile, probs = 0.975)
)

# Plot with smooth CI region
ggplot(predictions, aes(x = z_BASE_VELOCITY)) +
  # Add filled area for uncertainty
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.2) +
  # Add mean line
  geom_line(aes(y = estimate), color = "blue", size = 0.8) +
  # Add original data points
  geom_point(data = data, aes(y = FALLER), alpha = 0.5, size = 1) +
  # Scale adjustments
  scale_y_continuous(limits = c(0, 1)) +
  # Labels and theme
  labs(x = "Standardized Walking Velocity",
       y = "Probability of Fall") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")
  )
