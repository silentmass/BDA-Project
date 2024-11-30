data <- readRDS("data/clinical_demog_clean.rds")

# data$FALLER <- factor(data$FALLER)

data <- data %>%
  mutate(AGE_GROUP = case_when(
    AGE < 75 ~ 1,
    AGE >= 75 & AGE <= 80 ~ 2,
    AGE > 80 ~ 3
  ))

# Log transform TMT and z-score scale and fit

# Log transform TMTs
data$log_TMT_A <- log(data$TMT_A)
data$log_TMT_B <- log(data$TMT_B)

# Set original predictor categories
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

# Scale all predictors including log-transformed TMT variables
data[paste0("z_", original_predictors)] <- scale(data[original_predictors]) # Scale also GENDER


# Get total counts by gender and faller status
summary_table <- table(data$GENDER, data$FALLER)

# Convert to data frame and add readable labels
summary_df <- as.data.frame(summary_table)
colnames(summary_df) <- c("Gender", "Faller", "Count")
summary_df$Gender <- factor(summary_df$Gender, 
                            levels = c(0, 1),
                            labels = c("Male", "Female"))
summary_df$Faller <- factor(summary_df$Faller,
                            levels = c(0, 1), 
                            labels = c("Non-faller", "Faller"))

mean(data$AGE)
sd(data$AGE)
min(data$AGE)
max(data$AGE)

# Calculate percentages
total <- sum(summary_df$Count)
summary_df$Percentage <- round(summary_df$Count/total * 100, 1)

# Print results
print("Total counts by gender and faller status:")
print(summary_df)

print(paste("\nTotal number of subjects:", total))