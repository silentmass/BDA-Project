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

if (!require(shinystan)) {
  install.packages("shinystan")
  library(shinystan)
}

if (!require(latex2exp)) {
  install.packages("latex2exp")
  library(latex2exp)
}
pckgs_installed <- installed.packages()[,"Package"]
if (! "cspplot" %in% pckgs_installed) {
  remotes::install_github("CogSciPrag/cspplot")
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
fit_prefix <- "faller_classification_top2"

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
fits <- readRDS(paste0("models/", fit_prefix, "/faller_classification_models.rds"))

# Plot pp_check ----
set.seed(SEED)

# Filter the dataframe for selected models
fit_name <- "speed_dep_hier_gender_v2"
fit <- fits[[fit_name]]

subtitle <- TeX(paste0("Hierarchical SPEED DEPRESSION grouped by GENDER v2"))

nd <- nuts_params(fit) |>
  subset(Parameter == "divergent__") |>
  pull(Value) |>
  sum()

paste("Divergent transitions:", nd)

pp_plot <- pp_check(fit, ndraws = 50) +
  labs(title = "Posterior Predictive Check", subtitle = subtitle,
       x = "Probability",
       y = "Density",
       caption = "Dark line shows observed data. Light blue lines show 50 posterior predictions."
       ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray30"),
    plot.caption = element_text(size = 10, color = "gray30"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "white", linewidth = 0)
  ) +
  coord_cartesian(ylim = c(0, 2.0)) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25))

pp_plot

ggsave(
  paste0(c(plots_path, paste0("pp_check_", fit_name, ".png")), collapse = "/"),
  pp_plot,
  width = 7,     # inches
  height = 5,     # inches
  dpi = 300,      # high resolution
  bg = "white"    # ensure white background
)


# Plot factors ----
set.seed(SEED)

# global color scheme from CSP
project_colors = cspplot::list_colors() |> pull(hex)

plot_factors <- fit |> 
  tidy_draws() |> 
  select(starts_with("b_")) |>
  pivot_longer(cols = everything()) |> 
  mutate(
    name = str_remove(name, "b_z_"),
    name = str_remove(name, "b_"),
    name = case_match(
      name,
      "GDS" ~ "Depression (GDS)",
      "BASE_VELOCITY" ~ "Walking Speed",
      "Intercept" ~ "Intercept"
    ),
    name = reorder(name, abs(value), median)
  ) |>
  ggplot(aes(x = value, y = name)) +
  tidybayes::stat_halfeye(
    fill = project_colors[1],
    alpha = 0.7,
    .width = c(0.5, 0.95)  # Specify probability mass to show
  ) +
  labs(
    title = "Factors Associated with Fall Risk",  # Adjusted TeX formatting
    subtitle = subtitle,
    x = "Coefficient Estimate",
    y = "",
    caption = "Note: Higher values indicate increased fall risk. GDS = Geriatric Depression Scale"
  ) +
  geom_vline(xintercept = 0, color = project_colors[2], 
             size = 0.5, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray30"),
    plot.caption = element_text(size = 10, color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "white", linewidth = 0)
  ) +
  scale_x_continuous(
    limits = c(-2.5, 2.5),  # Adjusted limits to avoid cutting off distributions
    breaks = seq(-2, 2, by = 1)
  )

plot_factors

ggsave(
  paste0(c(plots_path, paste0("factors_", fit_name, ".png")), collapse = "/"),
  plot_factors,
  width = 7,     # inches
  height = 5,     # inches
  dpi = 300,      # high resolution
  bg = "white"    # ensure white background
)
