######## Run full ----

# Load librariers ----
library(ggplot2)
library(tidyr)
library(dplyr)

# load data ----
loos_normal_prior_zero_one <- read.csv("results/faller_classification_normal-prior_0-1/loo-comparison.csv")
loos_normal_prior_zero_five <- read.csv("results/faller_classification_normal-prior_0-5/loo-comparison.csv")
loos_studentt_prior_three <- read.csv("results/faller_classification_student-t-prior_3/loo-comparison.csv")
loos_studentt_prior_five <- read.csv("results/faller_classification_student-t-prior_5/loo-comparison.csv")

# Plot ----

# Add identifier column to each dataset
loos_normal_prior_zero_one$prior <- "Normal_0-1"
loos_normal_prior_zero_five$prior <- "Normal_0-5"
loos_studentt_prior_three$prior <- "Student-t_3"
loos_studentt_prior_five$prior <- "Student-t_5"

# Combine datasets
combined_loos <- rbind(
  loos_normal_prior_zero_one,
  loos_normal_prior_zero_five,
  loos_studentt_prior_three,
  loos_studentt_prior_five
)

# Clean model names
combined_loos$model <- gsub("c-", "", combined_loos$model)
combined_loos$model <- gsub("_", " ", combined_loos$model)

# Add this before the ggplot command:
# First set the factor levels for all models
model_order <- c(
  "PHYSICAL SPEED COGNITIVE DEPRESSION",
  "PHYSICAL SPEED",
  "PHYSICAL",
  "SPEED",
  "COGNITIVE",
  "spline-COGNITIVE",
  "DEPRESSION",
  "hierarchical-AGE GROUP-PHYSICAL SPEED COGNITIVE DEPRESSION",
  "hierarchical-GENDER-PHYSICAL SPEED COGNITIVE DEPRESSION",
  "hierarchical-AGE GROUP-SPEED",
  "hierarchical-GENDER-SPEED DEPRESSION",
  "hierarchical-AGE GROUP-spline-COGNITIVE"
)
combined_loos$model <- factor(combined_loos$model, levels = rev(model_order))

# Define selected models
selected_models <- c(
  "PHYSICAL SPEED COGNITIVE DEPRESSION",
  "SPEED",
  "DEPRESSION",
  "hierarchical-AGE GROUP-PHYSICAL SPEED COGNITIVE DEPRESSION",
  "hierarchical-GENDER-PHYSICAL SPEED COGNITIVE DEPRESSION",
  "hierarchical-AGE GROUP-SPEED",
  "hierarchical-GENDER-SPEED DEPRESSION"
)

# Filter the dataframe for selected models
plot_data <- combined_loos[combined_loos$model %in% selected_models, ]

# Then your existing ggplot code
plot_comp <- ggplot(plot_data, aes(x = elpd, y = model, color = prior)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = elpd - se, xmax = elpd + se),
                 position = position_dodge(width = 0.5),
                 height = 0.2) +
  theme_minimal() +
  labs(x = "ELPD", y = "Model", title = "LOO Comparison: Normal vs Student-t Priors") +
  theme(axis.text.y = element_text(size = 8)) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", linewidth = 0)
  )

ggsave(paste0(c("plots", paste0("large_models_sensitivity_comparison.png")), collapse = "/"), plot_comp, width = 10, height = 6)
