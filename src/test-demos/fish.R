n_draw <- 10000
n_data <- 16
original_subscribers <- 6

# draw from uniform prior
prior_rate <- runif(n_draw, 0, 1)

gen_model <- function(rate, n_data) {
  subscribers <- rbinom(1, size = n_data, prob = rate)
  return(subscribers)
}

subscribers <- rep(NA, n_draw)
for (i in 1:n_draw) {
  subscribers[i] <- gen_model(prior_rate[i], n_data)
}

post_rate <- prior_rate[subscribers == original_subscribers]

length(post_rate)

hist(post_rate, xlim = c(0, 1))