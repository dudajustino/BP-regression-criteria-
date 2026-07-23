rm(list=ls())

setwd("C:/Users/euedu/Downloads/git hub press/source codes")

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R")
source("criteria.R")

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

############### Corn dataset ##########################

#The data provide measurements of diameter (in inches at 4.5 feet above ground), height (feet),
#and volume (cubic feet) of wood in 31 black cherry trees felled in the Allegheny National Forest in Pennsylvania.
#Ryan et al. (1976, p. 329).

corn <- read.table("milho.dat", header = TRUE)

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
