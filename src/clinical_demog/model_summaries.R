######## Run full ----

# Init ----

#### Load packages and set theme ----

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

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
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

if (!require(projpred)) {
  install.packages("projpred")
  library(projpred)
}

cmdstan_installed <- function(){
  res <- try(out <- cmdstanr::cmdstan_path(), silent = TRUE)
  !inherits(res, "try-error")
}

if(!cmdstan_installed()){
  install_cmdstan()
}

ggplot2::theme_set(ggplot2::theme_minimal())

#### Set SEED ----
fit_prefix <- "faller_classification_normal-prior_0-1"

results_path <- paste0(c("results", fit_prefix), collapse = "/")
models_path <- paste0(c("models", fit_prefix), collapse = "/")
plots_path <- paste0(c("plots", fit_prefix), collapse = "/")
mcmc_path <- paste0(c("plots", fit_prefix, "MCMC"), collapse = "/")
mcmc_areas_path <- paste0(c("plots", fit_prefix, "MCMC", "areas"), collapse = "/")
mcmc_traces_path <- paste0(c("plots", fit_prefix, "MCMC", "traces"), collapse = "/")
mcmc_intervals_path <- paste0(c("plots", fit_prefix, "MCMC", "intervals"), collapse = "/")

SEED = 2024

set.seed(SEED)


#### Source helper functions ----
utils_path <- paste0(c("src", "clinical_demog", "utils"), collapse = "/")

source(paste0(c(utils_path, "analysis_plotting_helper_functions.R"), collapse = "/"))
source(paste0(c(utils_path, "analysis_helper_functions.R"), collapse = "/"))

# Load models
fits <- readRDS("models/faller_classification_top2/faller_classification_models.rds")

library(brms)
library(tidyverse)

get_prior_info <- function(fit, model_name) {
  prior_df <- as.data.frame(prior_summary(fit))
  prior_df$model <- model_name
  return(prior_df)
}

# Apply the function to your list of fits and combine results
prior_df <- do.call(rbind, 
                    Map(get_prior_info, 
                        fits, 
                        names(fits)))

# write.csv(prior_df, paste0(c(models_path, "model_priors.csv"), collapse = "/"), row.names = FALSE)

# Extract summaries ----

extract_speed_model_summaries <- function(fits) {
  # Helper function to safely extract columns that might not exist
  safe_extract <- function(matrix, col_name) {
    if(col_name %in% colnames(matrix)) {
      return(matrix[, col_name])
    } else {
      return(rep(NA, nrow(matrix)))
    }
  }
  
  # Extract fixed effects from all models
  fixed_effects <- map_dfr(names(fits), function(model_name) {
    model <- fits[[model_name]]
    fixed <- fixef(model)
    
    data.frame(
      model = model_name,
      parameter = rownames(fixed),
      estimate = safe_extract(fixed, "Estimate"),
      error = safe_extract(fixed, "Est.Error"),
      lower_ci = safe_extract(fixed, "Q2.5"),
      upper_ci = safe_extract(fixed, "Q97.5"),
      bulk_ess = safe_extract(fixed, "Bulk_ESS"),
      tail_ess = safe_extract(fixed, "Tail_ESS"),
      rhat = safe_extract(fixed, "Rhat"),
      stringsAsFactors = FALSE
    )
  })
  
  # Extract random effects if they exist
  random_effects <- map_dfr(names(fits), function(model_name) {
    model <- fits[[model_name]]
    
    tryCatch({
      if(!is.null(ranef(model))) {
        rand_summary <- VarCorr(model)
        rand_df <- map_dfr(names(rand_summary), function(group) {
          summary_mat <- rand_summary[[group]]
          data.frame(
            model = model_name,
            group = group,
            parameter = rownames(summary_mat),
            estimate = safe_extract(summary_mat, "Estimate"),
            error = safe_extract(summary_mat, "Est.Error"),
            lower_ci = safe_extract(summary_mat, "Q2.5"),
            upper_ci = safe_extract(summary_mat, "Q97.5"),
            bulk_ess = safe_extract(summary_mat, "Bulk_ESS"),
            tail_ess = safe_extract(summary_mat, "Tail_ESS"),
            rhat = safe_extract(summary_mat, "Rhat"),
            stringsAsFactors = FALSE
          )
        })
        return(rand_df)
      }
      return(NULL)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Extract model fit metrics
  fit_metrics <- map_dfr(names(fits), function(model_name) {
    model <- fits[[model_name]]
    
    # Get Bayesian R² with summary statistics
    r2_summary <- bayes_R2(model)
    
    # Get model summary
    model_summary <- summary(model)
    
    data.frame(
      model = model_name,
      r2_estimate = mean(r2_summary),
      r2_error = sd(r2_summary),
      r2_lower = quantile(r2_summary, probs = 0.025),
      r2_upper = quantile(r2_summary, probs = 0.975),
      looic = tryCatch(loo(model)$estimates["looic", "Estimate"], error = function(e) NA),
      waic = tryCatch(waic(model)$estimates["waic", "Estimate"], error = function(e) NA),
      stringsAsFactors = FALSE
    )
  })
  
  # Return all summaries in a list
  list(
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    fit_metrics = fit_metrics
  )
}

summaries <- extract_speed_model_summaries(fits)

# Check the results
head(summaries$fixed_effects)
head(summaries$random_effects)
head(summaries$fit_metrics)