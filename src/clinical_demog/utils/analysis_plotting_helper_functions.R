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