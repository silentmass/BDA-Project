# Load packages and set theme ----

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

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(pROC)) {
  install.packages("pROC")
  library(pROC)
}

if (!require(patchwork)) {
  install.packages("patchwork")
  library(patchwork)
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

# BMRS fit models and view summaries and plots ----

#otin nyt malliin parametreja, joilla näytti plotin perusteella olevan vaikutusta
fall_pooled_formula <- bf(YEAR_FALL ~ 1 +
                            MMSE +
                            TUG +
                            EFI_EXEC_FUNC_INDEX +
                            GCS_NEUROTRAX,
                          family = "gaussian")

# Show default priors
get_prior(fall_pooled_formula, data = data)

fall_pooled_fit <- brm(
  formula = fall_pooled_formula,
  data = data
)

summary(fall_pooled_fit)

# Create default density plot
pp_check(fall_pooled_fit)


# Classify subjects ----
#BERG-ABC-TUG ----
fall_class_formula <- bf(FALLER ~ 1 +
                           BERG +
                         ABC_TOTAL_PERCENT +
                         TUG,
                         family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  # Intercept: weakly informative prior
  prior(student_t(3, 0, 2.5), class = "Intercept"),
 
  
  # BERG: negative association with falls (higher score = better balance)
  prior(normal(-0.5, 0.5), class = "b", coef = "BERG"),
  prior(lognormal(0, 1), class = "b", coef = "ABC_TOTAL_PERCENT"),
  prior(lognormal(0, 1), class = "b", coef = "TUG")
  
  # Velocity measures: both directions possible but likely small effect
 # prior(normal(0, 0.5), class = "b", coef = "BASE_VELOCITY"),
#  prior(normal(0, 0.5), class = "b", coef = "S3_VELOCITY")
)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = data,
  family = bernoulli(),
  prior = fall_class_priors
)

# Posterior predictive check
pp_check(fall_class_fit)
#----
#BERG only ----
fall_class_formula <- bf(FALLER ~ 1 +
                           BERG,
                           
                         family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  # Intercept: weakly informative prior
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  
  
  # BERG: negative association with falls (higher score = better balance)
  prior(normal(-1, 1), class = "b", coef = "BERG")

  
  # Velocity measures: both directions possible but likely small effect
  # prior(normal(0, 0.5), class = "b", coef = "BASE_VELOCITY"),
  #  prior(normal(0, 0.5), class = "b", coef = "S3_VELOCITY")
)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = data,
  family = bernoulli(),
  prior = fall_class_priors
)

# Posterior predictive check
pp_check(fall_class_fit)
#----
#depression----
fall_class_formula <- bf(FALLER ~ 1 +
                           GDS + GCS_NEUROTRAX,
                         
                         family = "bernoulli")

get_prior(fall_class_formula, data = data)

fall_class_priors <- c(
  # Intercept: weakly informative prior
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  
  
  # BERG: negative association with falls (higher score = better balance)
  prior(normal(0, 1), class = "b", coef = "GDS"),
  prior(normal(0, 1), class = "b", coef = "GCS_NEUROTRAX")
  
  
  # Velocity measures: both directions possible but likely small effect
  # prior(normal(0, 0.5), class = "b", coef = "BASE_VELOCITY"),
  #  prior(normal(0, 0.5), class = "b", coef = "S3_VELOCITY")
)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = data,
  family = bernoulli(),
  prior = fall_class_priors
)

# Posterior predictive check
pp_check(fall_class_fit)
#----
summary(fall_class_fit)
# Get predicted probabilities
predictions <- posterior_predict(fall_class_fit)
pred_probs <- colMeans(predictions)

# ROC curve
roc_obj <- roc(data$FALLER, pred_probs)
roc_plot <- ggroc(roc_obj) +
  labs(title = "ROC Curve for Fall Classification",
       subtitle = paste("AUC =", round(auc(roc_obj), 3))) +
  theme_minimal()

# Confusion matrix at 0.5 threshold
pred_class <- ifelse(pred_probs > 0.5, 1, 0)
conf_mat <- table(Predicted = pred_class, Actual = data$FALLER)

# Metrics
sensitivity <- conf_mat[2,2] / sum(conf_mat[,2])
specificity <- conf_mat[1,1] / sum(conf_mat[,1])
accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

# Plot results
print(roc_plot)
print("Confusion Matrix:")
print(conf_mat)
print(paste("Sensitivity:", round(sensitivity, 3)))
print(paste("Specificity:", round(specificity, 3)))
print(paste("Accuracy:", round(accuracy, 3)))

# Generate diagnostic plots
# plot_mcmc_diagnostics(fall_pooled_fit)

# Check convergence metrics
# check_convergence(fall_pooled_fit)

# Plot data variables ----

plot(data$AGE, data$FALLER)
plot(data$EFI_EXEC_FUNC_INDEX, data$FALLER)
plot(data$GCS_NEUROTRAX, data$FALLER)
plot(data$SIX_MONTHS_FALL, data$FALLER)
plot(data$YEAR_FALL, data$FALLER)
plot(data$ABC_TOTAL_PERCENT, data$YEAR_FALL)
plot(data$AGE, data$YEAR_FALL)
plot(data$EFI_EXEC_FUNC_INDEX, data$YEAR_FALL)
plot(data$GCS_NEUROTRAX, data$YEAR_FALL)
plot(data$SF36, data$YEAR_FALL)
plot(data$DGI, data$YEAR_FALL)
plot(data$BERG, data$YEAR_FALL)
plot(data$MMSE, data$YEAR_FALL)
plot(data$FSST, data$FALLER)
plot(data$TUG, data$YEAR_FALL)

## All variables: 


fall_class_formula <- bf(FALLER ~ 1 +GDS +
                           AGE + GCS_NEUROTRAX + EFI_EXEC_FUNC_INDEX + GENDER + 
                           ABC_TOTAL_PERCENT+ SF36 + MMSE + MOCA + FAB + 
                           TUG + FSST + BERG + DGI + TMT_A+ TMT_B,
                         
                         family = "bernoulli")

get_prior(fall_class_formula, data = data)


fall_class_priors <- c(
  # Intercept: weakly informative prior
  prior(normal(0,10), class = "Intercept"),
  
  
  # BERG: negative association with falls (higher score = better balance)
  prior(normal(0, 10), class = "b", coef = "GDS"),
  prior(normal(0, 10), class = "b", coef = "GCS_NEUROTRAX"),
  prior(normal(0, 10), class = "b", coef = "AGE"),
  prior(normal(0, 10), class = "b", coef = "EFI_EXEC_FUNC_INDEX"),
  prior(normal(0, 10), class = "b", coef = "GENDER"),
  prior(normal(0, 10), class = "b", coef = "ABC_TOTAL_PERCENT"),
  prior(normal(0, 10), class = "b", coef = "SF36"),
  prior(normal(0, 10), class = "b", coef = "MMSE"),
  prior(normal(0, 10), class = "b", coef = "MOCA"),
  prior(normal(0, 10), class = "b", coef = "FAB"),
  prior(normal(0, 10), class = "b", coef = "TUG"),
  prior(normal(0, 10), class = "b", coef = "FSST"),
  prior(normal(0, 10), class = "b", coef = "BERG"),
  prior(normal(0, 10), class = "b", coef = "DGI"),
  prior(normal(0, 10), class = "b", coef = "TMT_A"),
  prior(normal(0, 10), class = "b", coef = "TMT_B")

)

fall_class_fit <- brm(
  formula = fall_class_formula,
  data = data,
  family = bernoulli(),
  prior = fall_class_priors
)

# Posterior predictive check
pp_check(fall_class_fit)
#----
summary(fall_class_fit)
# Get predicted probabilities