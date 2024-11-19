cStick = read.csv("cStick.csv", sep = ',')
par(mfrow = c(2, 3))
plot(data$Pressure, data$Decision)
plot(data$Sugar.level, data$Decision)
plot(data$HRV, data$Decision)
plot(data$SpO2, data$Decision)
plot(data$Distance, data$Decision)
plot(data$Accelerometer, data$Decision)


if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require(tidybayes)) {
  install.packages("tidybayes")
  library(tidybayes)
}

if (!require(brms)) {
  install.packages("brms")
  library(brms)
}

if (!require(metadat)) {
  install.packages("metadat")
  library(metadat)
}

if(!require(cmdstanr)){
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  library(cmdstanr)
}
cmdstan_installed <- function(){
  res <- try(out <- cmdstanr::cmdstan_path(), silent = TRUE)
  !inherits(res, "try-error")
}
if(!cmdstan_installed()){
  install_cmdstan()
}
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cStick = read.csv("data/cStick.csv", sep = ',')
fall_formula <- bf(
  Decision ~ 1 + Distance + HRV + SpO2,
  family = "gaussian",
  center = FALSE
)
get_prior(fall_formula, data = cStick)  

fall_fit <- brm(
  formula = fall_formula,
  data = cStick
)
summary(fall_fit)
pp_check(fall_fit)
