
####### Run pre ----

#### Load models and libraries ----

# prediction_script.R
library(brms)
library(tidyverse)
library(pROC)
library(ggdist)

# Set SEED ----
SEED = 2024
set.seed(SEED)

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
  BASE_VELOCITY = seq(from = min(data$BASE_VELOCITY), 
                      to = max(data$BASE_VELOCITY), 
                      length.out = 100),
  AGE = c(70, 80)
)

# Scale predictors
new_data$z_BASE_VELOCITY <- scale(new_data$BASE_VELOCITY)[,1]
new_data$z_AGE <- scale(new_data$AGE)[,1]

# Get both types of predictions
# Get predictions for probability curve
preds_epred <- posterior_epred(fit, newdata = new_data)
# Get one set of binary predictions (not averaging them)
preds_predict <- posterior_predict(fit, newdata = new_data)
# Take just one draw from the posterior for binary predictions
binary_preds <- preds_predict[1,]

# Calculate proportion of 1s for each point (uncertainty)
binary_proportions <- colMeans(preds_predict)
# Calculate confidence intervals
binary_lower <- apply(preds_predict, 2, quantile, probs = 0.025)
binary_upper <- apply(preds_predict, 2, quantile, probs = 0.975)

predictions <- data.frame(
  new_data,
  # Mean probabilities (epred)
  estimate_epred = colMeans(preds_epred),
  lower_epred = apply(preds_epred, 2, quantile, probs = 0.025),
  upper_epred = apply(preds_epred, 2, quantile, probs = 0.975),
  
  estimate_predict = binary_proportions,
  lower_predict = binary_lower,
  upper_predict = binary_upper,
  # Binary predictions (single draw)
  pred_binary = binary_preds
)

# Plot
ggplot(predictions, aes(x = BASE_VELOCITY)) +
  geom_ribbon(aes(ymin = lower_epred, ymax = upper_epred), alpha = 0.2) +
  geom_line(aes(y = estimate_epred)) +
  geom_line(aes(y = estimate_predict)) +
  # Now these points will be either 0 or 1
  # geom_point(aes(x = BASE_VELOCITY, y = pred_binary), alpha = 0.5) +
  geom_dots(data = predictions, aes(y = pred_binary, x = BASE_VELOCITY, 
                                 side = ifelse(pred_binary, "bottom", "top")), 
            pch = 19, color = "grey20", scale = 0.05) +
  facet_wrap(~AGE) +
  ylim(0, 1) +
  labs(x = "Walking Velocity (m/s)",
       y = "Probability of being faller") +
  theme_minimal()

# ----

# Create prediction grid
raw_new_data <- expand.grid(
  BASE_VELOCITY = seq(from = min(data$BASE_VELOCITY), 
                    to = max(data$BASE_VELOCITY), 
                    length.out = 100),
  AGE = c(70, 80)
  # z_BASE_VELOCITY = seq(from = min(data$z_BASE_VELOCITY), 
  #                       to = max(data$z_BASE_VELOCITY), 
  #                       length.out = 100)
)

scaled_new_data <- raw_new_data
if ("TMT_A" %in% names(raw_new_data)) {
  scaled_new_data$TMT_A < log(raw_new_data$TMT_A)
}
scaled_new_data$BASE_VELOCITY <- scale(raw_new_data$BASE_VELOCITY)[,1]
# scaled_new_data[] <- lapply(scaled_new_data, scale)  # Scale each column
names(scaled_new_data) <- paste0("z_", names(scaled_new_data))

# Get posterior predictions of probabilities
preds <- posterior_predict(fit, newdata = scaled_new_data)

# Calculate summary statistics
predictions <- data.frame(
  raw_new_data,
  # Mean prediction (probability)
  estimate = colSums(preds)/dim(preds)[1],
  # FALLER = round(colSums(preds)/dim(preds)[1]),
  # Uncertainty intervals
  lower = apply(preds, 2, quantile, probs = 0.025),
  upper = apply(preds, 2, quantile, probs = 0.975)
)

# Check original data distribution
ggplot(data, aes(x = BASE_VELOCITY, y = FALLER)) +
  geom_point(alpha = 0.5) +
  labs(x = "Walking Velocity (m/s)", y = "Faller Status") +
  theme_minimal()

# Plot with smooth CI region
ggplot(predictions, aes(x = BASE_VELOCITY)) +
  # Add filled area for uncertainty
  # geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.2) +
  # Add mean line
  geom_line(aes(y = estimate), color = "blue", size = 0.8) +
  # Add original data points
  geom_point(data = predictions, aes(y = estimate), alpha = 0.5, size = 1) +
  facet_wrap(~factor(AGE)) +
  # Scale adjustments
  scale_y_continuous(limits = c(0, 1)) +
  # Labels and theme
  labs(x = "Walking Velocity (m/s)",
       y = "Probability of being faller") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white")
  )
