library(ggplot2)
theme_set(theme_minimal())

n_girls <- 437
n_boys <- 543


# sequence of thetas
df1 <- data.frame(theta = seq(0.375, 0.525, 0.001))

a <- n_girls + 1 # Girls
b <- n_boys + 1 # Boys


# get posterior
df1$p <- dbeta(df1$theta, a, b)

df2 <- data.frame(theta = seq(qbeta(0.025, a, b), qbeta(0.975, a, b), length.out=100))

# get posterior from quantiles thetas
df2$p <- dbeta(df2$theta, a, b)