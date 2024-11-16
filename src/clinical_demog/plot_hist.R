library(ggplot2)
library(dplyr)  # for handling tibbles

# Read the data
data <- readRDS("data/clinical_demog_clean.rds")

# Get numeric columns
numeric_cols <- names(data)[sapply(data, is.numeric)]
numeric_cols <- numeric_cols[numeric_cols != "FALLER"]

# Create and save histogram for each numeric column
for(col in numeric_cols) {
  # Get values using pull() for tibbles
  stats_0 <- data %>% filter(FALLER == 0) %>% pull(col)
  stats_1 <- data %>% filter(FALLER == 1) %>% pull(col)
  
  stats_text <- paste(
    "Non-fallers (0):",
    "\nCount:", length(stats_0),
    "\nMean ± SD:", round(mean(stats_0, na.rm = TRUE), 2), "±", round(sd(stats_0, na.rm = TRUE), 2),
    "\nMedian [IQR]:", round(median(stats_0, na.rm = TRUE), 2), 
    "[", round(quantile(stats_0, 0.25, na.rm = TRUE), 2), "-", round(quantile(stats_0, 0.75, na.rm = TRUE), 2), "]",
    "\n\nFallers (1):",
    "\nCount:", length(stats_1),
    "\nMean ± SD:", round(mean(stats_1, na.rm = TRUE), 2), "±", round(sd(stats_1, na.rm = TRUE), 2),
    "\nMedian [IQR]:", round(median(stats_1, na.rm = TRUE), 2),
    "[", round(quantile(stats_1, 0.25, na.rm = TRUE), 2), "-", round(quantile(stats_1, 0.75, na.rm = TRUE), 2), "]"
  )
  
  p <- ggplot(data, aes(x = .data[[col]], fill = factor(FALLER))) +
    geom_histogram(position = "dodge", alpha = 0.8, bins = 30) +
    labs(title = paste("Histogram of", col),
         x = col,
         fill = "FALLER") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white")
    ) +
    scale_fill_manual(values = c("0" = "#2166AC", "1" = "#B2182B")) +
    annotate("text", 
             x = max(data[[col]], na.rm = TRUE),  # Changed to max for right positioning
             y = Inf,
             label = stats_text,
             hjust = 1,  # Changed to 1 for right alignment
             vjust = 1,
             size = 3.5,
             fontface = 2)
  
  # Save plot
  ggsave(
    filename = paste0("plots/histograms_clean/", col, "_clean_histogram.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  print(p)
}