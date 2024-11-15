load("/Users/juha/BDA_course_Aalto/rpackage/data/algae.rda")
library(ggplot2)
library(glue)

get_posterior_mean <- function(alpha_posterior, beta_posterior) {
  mean_E <- alpha_posterior / (alpha_posterior + beta_posterior)
  return(mean_E)
}

get_prior_and_posterior <- function(theta_values, prior_a, prior_b) {
  
  alpha_posterior <- prior_a+y
  beta_posterior <- prior_b+n-y
  
  prior <- dbeta(theta_values, prior_b, prior_b)
  posterior <- dbeta(theta_values, alpha_posterior, beta_posterior)
  mean_E <- get_posterior_mean(alpha_posterior, beta_posterior)
  
  return(
    list(
      prior = prior,
      posterior = posterior,
      mean_E = mean_E,
      alpha_posterior = alpha_posterior,
      beta_posterior = beta_posterior
    )
  )
}

get_credibel_interval <- function(alpha_posterior, beta_posterior, lo_bound = 0.05, hi_bound = 0.95) {
  # For a 90% credible interval
  lower <- round(qbeta(lo_bound, alpha_posterior, beta_posterior), 2)
  upper <- round(qbeta(hi_bound, alpha_posterior, beta_posterior), 2)
  
  return(list(
    lower = lower,
    upper = upper,
    lo_bound = lo_bound,
    hi_bound = hi_bound
  ))
}

get_interval_df <- function(alpha_posterior,
                            beta_posterior,
                            lo_bound = 0.025,
                            hi_bound = 0.975) {
  # seq creates evenly spaced values from 2.5% quantile
  # to 97.5% quantile (i.e., 95% central interval)
  # qbeta computes the value for a given quantile given parameters a and b
  df <- data.frame(theta = seq(
    qbeta(lo_bound, alpha_posterior, beta_posterior),
    qbeta(hi_bound, alpha_posterior, beta_posterior),
    length.out = 100
  ))
  # compute the posterior density
  df$posterior <- dbeta(df$theta, alpha_posterior, beta_posterior)
  return(df)
}

get_beta_posterior_less_or_equal_to_thres <- function(thres, alpha_posterior, beta_posterior) {
  return(pbeta(thres, alpha_posterior, beta_posterior))
}

plot_geom_jitter_binomial_data <- function(df, x, y) {
  ggplot(df, aes(x = x, y = y)) +
    geom_jitter(height = 0.05, alpha = 0.5) +
    scale_y_continuous(limits = c(-0.2, 1.2))
}


df_algae <- data.frame(lake = seq(0, 1, length.out = length(algae)), y = algae)

y <- sum(algae)
n <- length(algae)

# historical detection rate
theta_zero <- 0.2

theta_values <- seq(0, 1, length.out = 1000)

data1 <- get_prior_and_posterior(theta_values = theta_values, prior_a = 1, prior_b = 1)

df1 <- data.frame(theta = theta_values, prior = data1$prior, posterior = data1$posterior)

credibel_intervals <- get_credibel_interval(alpha_posterior, beta_posterior)

df_interval <- get_interval_df(
  alpha_posterior = alpha_posterior, 
  beta_posterior = beta_posterior
)

get_beta_posterior_less_or_equal_to_thres(
  theta_zero, alpha_posterior, beta_posterior
)


ggplot() + 
  geom_line(mapping = aes(theta, prior, color = "Prior"), data = df1) + 
  geom_line(data = df1, mapping = aes(theta, posterior, color = "Posterior")) +
  geom_area(data = df_interval, mapping = aes(theta, posterior, color = "Interval", fill = '1')) + 
  geom_vline(xintercept = theta_zero, linetype='dotted') + 
  labs(
    title = glue("Prior Beta({prior_a},{prior_b}) -> Posterior is Beta({alpha_posterior},{beta_posterior})"),
    x = "Theta",
    y = "Density",
    color = "Distribution"  # Legend title
  ) + 
  scale_y_continuous(expand = c(0, 0.1), breaks = NULL) +
  scale_fill_manual(values = 'lightblue', labels = '95% posterior interval') +
  theme(legend.position = 'bottom', legend.title = element_blank())

# Jittered points



plot_geom_jitter_binomial_data(df = df_algae, x = "lake", y = "y")
