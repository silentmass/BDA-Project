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

get_mcmc_theme <- function() {
  mcmc_fontsize <- 12
  mcmc_theme <- theme_minimal() +
    theme(
      text = element_text(size = mcmc_fontsize),
      axis.title = element_text(size = mcmc_fontsize),
      axis.text = element_text(size = mcmc_fontsize),
      strip.text = element_text(size = mcmc_fontsize, face = "bold"),
      legend.text = element_text(size = mcmc_fontsize),
      legend.title = element_text(size = mcmc_fontsize),
      panel.border = element_rect(color = "grey90", 
                                  fill = NA),
      plot.background = element_rect(fill = "white", linewidth = 0),
      panel.spacing.x = unit(0.7, "cm"),
      panel.spacing.y = unit(0.5, "cm"),
      aspect.ratio = 0.7,
      strip.text.y = element_text(angle = 0, 
                                  face = "bold", 
                                  size = mcmc_fontsize,
                                  margin = margin(r = 10)),
      strip.text.x = element_text(margin = margin(b = 10))
    )
  return(mcmc_theme)
}

plot_mcmc <- function(fit_object, fit_name, predictors) {
  posterior <- as_draws_df(fit_object)
  mcmc_theme <- get_mcmc_theme()
  
  mcmc_areas(
    posterior,
    pars = paste0("b_", predictors)
  ) + mcmc_theme
  ggsave(paste0(c("plots", "predictors", paste0(paste0(c("mcmc-areas", fit_name), collapse = "_"), ".png")), collapse = "/"))
  
  mcmc_trace(
    posterior,
    pars = paste0("b_", predictors)
  ) + mcmc_theme
  ggsave(paste0(c("plots", "predictors", paste0(paste0(c("mcmc-traces", fit_name), collapse = "_"), ".png")), collapse = "/"))
  
  mcmc_intervals(
    posterior,
    pars = paste0("b_", predictors),
  ) + mcmc_theme
  ggsave(paste0(c("plots", "predictors", paste0(paste0(c("mcmc-intervals", fit_name), collapse = "_"), ".png")), collapse = "/"))
}