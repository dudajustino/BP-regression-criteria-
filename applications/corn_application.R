rm(list=ls())

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R")
source("criteria.R")

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

############### Corn dataset ##########################

# Dataset available in Griffiths et al. (1993) on corn yield (pounds/acre)
# as a function of combinations of nitrogen and phosphate fertilizer levels.
#
# Griffiths, W. E., Hill, R. C., & Judge, G. G. (1993).
# Learning and Practicing Econometrics. Wiley.

#X1: nitrogen level
#X2: phosphate level
#Y : corn yield (pounds/acre)

corn <- read.table("milho.dat", header = TRUE)
head(corn)

# candidate models

fit <- gamlss(prod~I(log(ni)),family = BP(mu.link = "log"), data = corn, trace=FALSE)
summary(fit)

fit1 <- gamlss(prod~I(log(ni))+I(log(fo)),family = BP(mu.link = "sqrt"), data = corn, trace=FALSE)
summary(fit1)

fit2 <- gamlss(prod~I(log(ni))+I(log(fo)),family = BP(mu.link = "log"), data = corn, trace=FALSE)
summary(fit2)

fit3 <- gamlss(prod~I(log(ni))+I(log(fo)),~fo,family = BP(mu.link = "log", sigma.link = "log"), data = corn, trace=FALSE)
summary(fit3)

xtable(rbind(criteria(fit),
             criteria(fit1),
             criteria(fit2),
             criteria(fit3)), digits = 3)

# fit3 - model selected by the predictive and goodness-of-fit measures
envelope.BP(fit3, k=100, type = "quantile", link=c("log","log"))
plot.BP(fit3, which = 1, type = "quantile", is_application=T, q1=-2,q2=2, pos2 = c(3))
plot.BP(fit3, which = 2, type = "quantile", is_application=T, q1=-2,q2=2, pos2 = c(3))

# LOOCV procedure
n <- nrow(corn)
y_obs <- corn$prod
y_pred <- numeric(n)

for(i in 1:n){

  train_data <- corn[-i, ]
  test_data  <- corn[i, , drop = FALSE]

  fit_cv <- gamlss(
    prod ~ I(log(ni))+I(log(fo)),
    sigma.formula = ~ fo,
    family = BP(mu.link = "log", sigma.link = "log"),
    trace = FALSE,
    data = train_data
  )

  # previsão da média mu para a observação deixada de fora
  mu_pred <- predict(fit_cv, newdata = test_data, what = "mu", type = "response")

  y_pred[i] <- mu_pred
}

erro <- y_obs - y_pred

mspe <- mean(erro^2)
mape <- mean(abs(erro / y_obs)) * 100
mae  <- mean(abs(erro))
rmse <- sqrt(mspe)

result <- list(
  y_obs = y_obs,
  y_pred = y_pred,
  erro = erro,
  mspe = mspe,
  rmse = rmse,
  mae = mae,
  mape = mape
)

result
