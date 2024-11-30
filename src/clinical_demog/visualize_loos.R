# load data ----

loos_normal_priors <- read.csv("results/faller_classification/loo-comparison.csv")
loos_studentt_priors <- read.csv("results/faller_classification_student-t-prior/loo-comparison.csv")


# Plot ----
library(ggplot2)
library(tidyr)
library(dplyr)

# Add identifier column to each dataset
loos_normal_priors$prior <- "Normal"
loos_studentt_priors$prior <- "Student-t"

# Combine datasets
combined_loos <- rbind(loos_normal_priors, loos_studentt_priors)

# Clean model names
combined_loos$model <- gsub("c-", "", combined_loos$model)
combined_loos$model <- gsub("_", " ", combined_loos$model)

# Create plot
ggplot(combined_loos, aes(x = elpd, y = model, color = prior)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = elpd - se, xmax = elpd + se),
                 position = position_dodge(width = 0.5),
                 height = 0.2) +
  theme_minimal() +
  labs(x = "ELPD", y = "Model", title = "LOO Comparison: Normal vs Student-t Priors") +
  theme(axis.text.y = element_text(size = 8))