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

generate_pars <- function(model_name, predictors) {
  # Check if model is hierarchical
  is_hierarchical <- grepl("^hierarchical-", model_name)
  hierarchical_group <- if(is_hierarchical) {
    sub("^hierarchical-(.+?)-.*$", "\\1", model_name)
  } else {
    NULL
  }
  
  # Check if model includes splines
  has_splines <- grepl("spline", model_name)
  
  if(has_splines) {
    pars <- c(
      "b_Intercept",
      paste0("bs_sz_", gsub("z_", "", predictors), "_1"),
      paste0("sds_sz_", gsub("z_", "", predictors), "_1")
    )
    if(is_hierarchical) {
      pars <- c(pars, paste0("sd_", hierarchical_group, "__", predictors))
    }
  } else {
    pars <- paste0("b_", predictors)
    if(is_hierarchical) {
      pars <- c(pars, paste0("sd_", hierarchical_group, "__", predictors))
    }
  }
  
  return(pars)
}

format_model_name <- function(name) {
  name <- gsub("hierarchical", "h", name)
  name <- gsub("AGE_GROUP", "AGE", name)
  
  categories <- list()
  is_spline <- FALSE
  hierarchical <- list()
  
  if (grepl("c-", name)) {
    categories <- strsplit(strsplit(name, "c-")[[1]][2], "_")[[1]]
  }
  
  if (grepl("^h-|^hierarchical-", name)) {
    # get initial part from the name
    hierarchical <- strsplit(name, "c-")[[1]]
    hierarchical <- strsplit(name, "-")[[1]]
    hierarchical <- hierarchical[grepl("AGE_GROUP|AGE|GENDER", hierarchical)]
  }
  
  if (grepl("spline", name)) {
    # get initial part from the name
    is_spline <- TRUE
  }
  
  abbreviated_name <- ""
  
  if (is_spline) {
    abbreviated_name <- "Spline | "
  }
  
  if (length(hierarchical)) {
    abbreviated_name <- paste0(abbreviated_name, "Hierarchical: ", paste0(hierarchical, collapse = ", "), " | ")
  }
  
  if (length(categories)) {
    abbreviated_name <- paste0(abbreviated_name, "Categories: ", paste0(categories, collapse = ", "))
  }
  
  return(abbreviated_name)
}

# Main plotting function that handles all three types
plot_mcmc <- function(fit, model_name, predictors, plot_type = "area", 
                      fontsize = 14, aspect_ratio = 1.0) {
  
  posterior <- as_draws_df(fit)
  pars <- generate_pars(model_name, predictors)
  title <- get_plot_title(model_name, plot_type)
  
  base_theme <- theme_minimal(base_size = fontsize) +
    theme(
      plot.title = element_text(hjust = 0.5, size = fontsize, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = fontsize, color = "gray40"),
      axis.text.y = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(size = fontsize),
      legend.text = element_text(size = fontsize),
      legend.title = element_text(size = fontsize),
      plot.background = element_rect(fill = "white", linewidth = 0)
    )
  
  if(plot_type == "area") {
    p <- mcmc_areas(posterior, pars = pars) +
      labs(
        title = title,
        x = "Standardized Coefficient",
        y = "Parameter",
        subtitle = "Shaded areas represent 95% credible intervals"
      ) +
      base_theme +
      theme(aspect.ratio = aspect_ratio)
    
  } else if(plot_type == "trace") {
    p <- mcmc_trace(posterior, pars = pars) +
      labs(
        title = "MCMC Chain Traces",
        x = "Iteration",
        y = "Parameter Value",
        subtitle = format_model_name(model_name)
      ) +
      base_theme +
      theme(
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      ) +
      facet_wrap(~parameter, ncol = 4)
    
  } else if(plot_type == "interval") {
    p <- mcmc_intervals(posterior, pars = pars) +
      labs(
        title = title,
        subtitle = "95% Credible Intervals",
        x = "Standardized Effect",
        y = "Parameter"
      ) +
      base_theme +
      theme(aspect.ratio = aspect_ratio)
  }
  
  p + geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
}

# Helper function for titles
get_plot_title <- function(model_name, plot_type) {
  abbr_name <- format_model_name(model_name)
  switch(plot_type,
         "area" = paste("Posterior Distributions of", abbr_name),
         "trace" = "MCMC Chain Traces",
         "interval" = paste("Posterior Intervals of", abbr_name)
  )
}