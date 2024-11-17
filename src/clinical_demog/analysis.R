install.packages("brms")
library(brms)

data <- readRDS("data/clinical_demog_clean.rds")
head(data)

#tässä on ensimmäinen malli, ei näytä toimivan, eli formulaa pitää muuttaa. 
#nyt siinä on vain yksi parametri (AGE), jolla ei summaryn perusteella ole vaikutusta
#periaatteessa olisi helppo tehdä kaksi eri mallia samoilla parametreilla: 
#toinen pooled ja toinen hierarchical
fall_pooled_formula <- bf(
  FALLER | trials(1) ~ 1 + AGE,
  family = binomial()
  )

fall_pooled_fit <- brm(
  formula = fall_pooled_formula,
  data = data
)

summary(fall_pooled_fit)