rm(list=ls())

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R")
source("criteria.R")

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

############### Land rent data dataset ##########################

#dataset available at Weisberg [28] on the average rent paid (in US dollars) per acre planted to alfalfa in 67
#counties in Minnesota in 1977
#Weisberg, S. (2014). Applied Linear Regression, 4th edition. Hoboken NJ: Wiley.

library("alr4")

#X1: average rent paid for land for other agricultural purposes
#X2: density of dairy cows in number per square mile
#X3: proportion of agricultural farmland used for pasture
#X4: 1 if liming is necessary for alfalfa cultivation; 0 otherwise
#Y: average rent paid per acre planted to alfalfa

data(landrent)
head(landrent)
landrent <- subset(landrent, X4 != 0)

# candidate models

fit <- gamlss(Y~I(log(X1)),family = BP(mu.link = "log"), data = landrent,trace=FALSE)
summary(fit)

fit1 <- gamlss(Y~I(log(X1))+I(log(X2)),family = BP(mu.link = "sqrt"), data = landrent,trace=FALSE)
summary(fit1)

fit2 <- gamlss(Y~I(log(X1))+I(log(X2)),family = BP(mu.link = "log"), data = landrent,trace=FALSE)
summary(fit2)

fit3 <- gamlss(Y~I(log(X1))+I(log(X2)),~X2,family = BP(mu.link = "log", sigma.link = "log"), data = landrent, trace=FALSE)
summary(fit3)

xtable(rbind(criteria(fit),
             criteria(fit1),
             criteria(fit2),
             criteria(fit3)), digits = 3)

# fit3 - model selected by the predictive and goodness-of-fit measures
envelope.BP(fit3, k=100, type = "quantile", link=c("log","log"))
plot.BP(fit3, which = 1, type = "quantile", is_application=T, q1=-2,q2=2, pos2 = c(3))
plot.BP(fit3, which = 2, type = "quantile", is_application=T, q1=-2,q2=2, pos2 = c(3))

# exclusion of the outlying observation (case 33)
landrent_33 <- landrent[-33,]
fit3_33 <- gamlss(Y~I(log(X1))+I(log(X2)),~X2,family = BP(mu.link = "log", sigma.link = "log"), data = landrent_33, trace=FALSE)
summary(fit3_33)

100 * ( c(fit3$mu.coefficients, fit3$sigma.coefficients) - c(fit3_33$mu.coefficients, fit3_33$sigma.coefficients) ) /
  (c(fit3$mu.coefficients, fit3$sigma.coefficients))

100 * (c(sqrt(diag(vcov(fit3)))) - c(sqrt(diag(vcov(fit3_33))))) / c(sqrt(diag(vcov(fit3))))

criteria(fit3)
criteria(fit3_33)

AIC(fit3, fit3_33)
BIC(fit3, fit3_33)

# LOOCV procedure
n <- nrow(landrent)
y_obs <- landrent$Y
y_pred <- numeric(n)

for(i in 1:n){

  train_data <- landrent[-i, ]
  test_data  <- landrent[i, , drop = FALSE]

  fit_cv <- gamlss(
    Y ~ I(log(X1)) + I(log(X2)),
    sigma.formula = ~ X2,
    family = BP(mu.link = "log", sigma.link = "log"),
    trace = FALSE,
    data = train_data
  )

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
