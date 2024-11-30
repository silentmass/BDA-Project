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

# Scaled predictor names
cols_list <- list(
  PHYSICAL = paste0("z_", original_cols_list[["PHYSICAL"]]),
  SPEED = paste0("z_", original_cols_list[["SPEED"]]),
  COGNITIVE = paste0("z_", original_cols_list[["COGNITIVE"]]),
  DEPRESSION = paste0("z_", original_cols_list[["DEPRESSION"]])
)

predictor_categories <- get_predictor_categories_from_cols_list(cols_list)

# Predict ----

# Prediction function
predict_faller_status <- function(model, new_data, threshold = 0.5) {
  predictions <- posterior_predict(model, newdata = new_data)
  pred_probs <- colMeans(predictions)
  pred_intervals <- apply(predictions, 2, quantile, probs = c(0.025, 0.975))
  
  data.frame(
    probability = pred_probs,
    lower_ci = pred_intervals[1,],
    upper_ci = pred_intervals[2,],
    predicted_class = ifelse(pred_probs > threshold, 1, 0)
  )
}

# Function to convert raw scores to z-scores
convert_to_zscores <- function(new_data, reference_data = data) {
  z_scores <- new_data %>%
    mutate(
      # Log transform TMTs first
      log_TMT_B = log(TMT_B),
      
      # Create z-scores for all variables
      z_AGE = scale(AGE, center = mean(reference_data$AGE), scale = sd(reference_data$AGE))[,1],
      z_GENDER = scale(GENDER, center = mean(reference_data$GENDER), scale = sd(reference_data$GENDER))[,1],
      z_DGI = scale(DGI, center = mean(reference_data$DGI), scale = sd(reference_data$DGI))[,1],
      z_TUG = scale(TUG, center = mean(reference_data$TUG), scale = sd(reference_data$TUG))[,1],
      z_FSST = scale(FSST, center = mean(reference_data$FSST), scale = sd(reference_data$FSST))[,1],
      z_BASE_VELOCITY = scale(BASE_VELOCITY, center = mean(reference_data$BASE_VELOCITY), scale = sd(reference_data$BASE_VELOCITY))[,1],
      z_log_TMT_B = scale(log_TMT_B, center = mean(reference_data$log_TMT_B), scale = sd(reference_data$log_TMT_B))[,1],
      z_GDS = scale(GDS, center = mean(reference_data$GDS), scale = sd(reference_data$GDS))[,1]
    ) %>%
    select(starts_with("z_"))
  
  return(z_scores)
}

# Example usage
new_raw_data <- data.frame(
  AGE = 78,
  GENDER = 1,
  DGI = 20,
  TUG = 12,
  FSST = 15,
  BASE_VELOCITY = 1.0,
  TMT_B = 120,
  GDS = 3
)

z_scored_data <- convert_to_zscores(new_raw_data)

# Get predictions
predictions <- predict_faller_status(fits[["c-PHYSICAL_SPEED"]], z_scored_data)

# Calculate ROC curve to find optimal threshold
roc_obj <- roc(data$FALLER, 
               posterior_predict(fits[["c-PHYSICAL_SPEED"]], newdata = data) %>% colMeans())
optimal_threshold <- coords(roc_obj, "best")["threshold"]

# new ----

new_person <- data.frame(
  AGE = 65,
  GENDER = 0,
  DGI = 20,
  TUG = 12,
  FSST = 15,
  BASE_VELOCITY = 1.0,
  TMT_B = 120,
  GDS = 3
)

z_scored_data <- convert_to_zscores(new_person)
result <- predict_faller_status(fits[["c-PHYSICAL_SPEED"]], z_scored_data, threshold = optimal_threshold)
print(result)

#### combinations ----

# Function to test different score combinations
test_combinations <- function(model, threshold = 0.520125, sample_size = 1000) {
  # Create smaller sequences
  test_cases <- expand.grid(
    AGE = seq(70, 80, by = 10),          # Reduced to 2 values
    GENDER = c(0, 1),
    DGI = seq(15, 24, by = 4.5),         # Reduced to 3 values
    TUG = seq(8, 16, by = 4),            # Reduced to 3 values
    FSST = seq(10, 20, by = 5),          # Reduced to 3 values
    BASE_VELOCITY = seq(0.8, 1.2, by = 0.2), # Reduced to 3 values
    TMT_B = c(110, 130),
    GDS = c(2, 4)
  )
  
  # If too many combinations, sample randomly
  if (nrow(test_cases) > sample_size) {
    test_cases <- test_cases[sample(nrow(test_cases), sample_size), ]
  }
  
  z_scored <- convert_to_zscores(test_cases)
  predictions <- predict_faller_status(model, z_scored, threshold)
  
  results <- cbind(test_cases, predictions)
  results[order(results$probability, decreasing = TRUE), ]
}

# Use function
results <- test_combinations(fits[["c-PHYSICAL_SPEED"]])
print(head(results))

# Basic visualization of probabilities across combinations
ggplot(results, aes(y = probability)) +
  geom_histogram() +
  labs(title = "Distribution of Fall Probabilities")

# Variable importance plot
ggplot(results, aes(x = reorder(FSST, probability), y = probability)) +
  geom_point() +
  facet_wrap(~GENDER) +
  labs(title = "Fall Probability by FSST Score and Gender", 
       x = "FSST Score", y = "Probability")

# Multiple variables heatmap
ggplot(results, aes(x = BASE_VELOCITY, y = FSST)) +
  geom_tile(aes(fill = probability)) +
  scale_fill_gradient(low = "blue", high = "red") +
  facet_grid(GENDER ~ AGE) +
  labs(title = "Fall Probability Heatmap")

# Boxplots for key metrics
results_long <- results %>%
  pivot_longer(c(DGI, TUG, FSST, BASE_VELOCITY), 
               names_to = "metric", values_to = "value")

# Add faceting to handle different scales appropriately
ggplot(results_long, aes(x = metric, y = value, fill = factor(predicted_class))) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free") +
  labs(title = "Raw Score Distributions by Predicted Class",
       x = "Metric", 
       y = "Score",
       fill = "Predicted Class") +
  scale_fill_discrete(labels = c("Non-faller", "Faller")) +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) # Add 20% padding above and below