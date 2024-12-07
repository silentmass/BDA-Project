library(ggplot2)

# Create data frame with x values
x <- seq(0, 3, length.out = 1000)
df <- data.frame(
  x = rep(x, 2),
  density = c(
    dexp(x, rate = 3),  # Exponential(3)
    dgamma(x, shape = 3, rate = 3)  # Gamma(3,3)
  ),
  Distribution = rep(c("Exponential(3)", "Gamma(3,3)"), each = length(x))
)

# Create plot
ggplot(df, aes(x = x, y = density, color = Distribution)) +
  geom_line(size = 1.2) +
  labs(
    title = "Comparison of Prior Distributions for SD Parameter",
    subtitle = "Exponential(3) vs Gamma(3,3)",
    x = "Standard Deviation",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("darkblue", "darkred")) +
  annotate("text", x = 1.5, y = 1.5, 
           label = "Exp(3): Mean = 0.33, Var = 0.11\nGamma(3,3): Mean = 1, Var = 0.33",
           size = 3.5)