source("src/clinical_demog/utils/analysis_helper_functions.R")

# Analysis plotting helper functions

# Function to create comprehensive MCMC diagnostics
plot_mcmc_diagnostics <- function(model) {
  # Extract posterior draws
  posterior_samples <- as.array(model)
  
  # Create individual plots
  p1 <- mcmc_trace(posterior_samples) +
    ggtitle("Trace Plots") +
    theme_minimal()
  
  p2 <- mcmc_dens_overlay(posterior_samples) +
    ggtitle("Parameter Densities") +
    theme_minimal()
  
  p3 <- mcmc_acf(posterior_samples, lags = 20) +
    ggtitle("Autocorrelation Plots") +
    theme_minimal()
  
  # Arrange plots
  p1 / p2 / p3
}

plot_pp_check <- function(fall_class_fit, selected_variables, variables_per_line = -1, ylim = c(0, 2.0)) {
  pp_check(fall_class_fit) +
    ggtitle(paste0("Posterior Predictive Check:\n", 
                   format_variables(selected_variables, variables_per_line = variables_per_line))) +
    ylim(ylim) +
    theme_minimal()
}