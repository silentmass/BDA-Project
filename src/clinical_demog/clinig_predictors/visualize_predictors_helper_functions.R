library(moments)
library(ggplot2)
library(patchwork)

compare_log_transformed <- function(df, col_name) {
  df[[paste0("log_", col_name)]] <- log(df[[col_name]])
  df_long <- df %>%
    pivot_longer(
      cols = all_of(c(col_name, paste0("log_", col_name))),
      names_to = "PREDICTOR",
      values_to = "PREDICTOR_VALUE"
    ) %>%
    select(ID, PREDICTOR, PREDICTOR_VALUE, FALLER)
  
  plot_data <- df_long %>%
    mutate(
      PREDICTOR = factor(PREDICTOR,
                         levels = c(col_name, paste0("log_", col_name))),
      FALLER = factor(FALLER)
    )
  
  # Common theme for all plots
  common_theme <- theme_minimal() +
    theme(
      text = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16)
    )
  
  # Calculate skewness and kurtosis for each predictor
  stats_df <- plot_data %>%
    group_by(PREDICTOR) %>%
    summarise(
      skew = round(skewness(PREDICTOR_VALUE), 2),
      kurt = round(kurtosis(PREDICTOR_VALUE), 2),
      x_pos = min(PREDICTOR_VALUE)
    )
  
  # Create histogram plot
  p1 <- ggplot(data = plot_data) +
    geom_histogram(mapping = aes(x = PREDICTOR_VALUE, 
                                 fill = FALLER,
                                 group = FALLER),
                   position = "identity", 
                   alpha = 0.3, 
                   bins = 30) +
    labs(y = str_to_upper("Count")) +
    facet_wrap(~PREDICTOR, scales = 'free_x',
               labeller = labeller(PREDICTOR = label_wrap_gen(width = 15)), 
               ncol=2) +
    geom_text(data = stats_df,
              aes(x = x_pos, 
                  y = Inf,
                  label = paste("Skewness:", skew, "\nKurtosis:", kurt)),
              hjust = 0, 
              vjust = 2,
              size = 5) +  # Adjusted text size for statistics
    common_theme
  
  # Create density plot
  p2 <- ggplot(data = plot_data) +
    geom_density(aes(x = PREDICTOR_VALUE,
                     color = FALLER,
                     fill = FALLER),
                 alpha = 0.3) +
    labs(y = str_to_upper("Density")) +
    facet_wrap(~PREDICTOR, scales = 'free',
               labeller = labeller(PREDICTOR = label_wrap_gen(width = 15)), 
               ncol=2) +
    common_theme
  
  # Create Q-Q plots
  p3 <- ggplot(data = plot_data) +
    geom_qq(aes(sample = PREDICTOR_VALUE, color = FALLER)) +
    geom_qq_line(aes(sample = PREDICTOR_VALUE, color = FALLER), 
                 linewidth = 0.8,
                 linetype = "dashed") +
    labs(y = str_to_upper("Sample Quantiles"), x = str_to_upper("Theoretical Quantiles")) +
    facet_wrap(~ PREDICTOR, scales = 'free_y',
               labeller = labeller(PREDICTOR = label_wrap_gen(width = 15)), 
               ncol=2) +
    common_theme
  
  # Combine plots using patchwork
  combined_plot <- p1 / p2 / p3 +
    plot_layout(guides = "collect", heights = c(1, 0.7, 0.7)) +
    plot_annotation(
      title = paste('Distribution Comparison:', col_name),
      theme = theme(
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(size = 16)
      )
    )
  
  return(combined_plot)
}