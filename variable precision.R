
# Clear all existing objects from the workspace
rm(list=ls())

setwd("C:/Users/euedu/Downloads/git hub press/source codes")

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R")
source("criteria.R")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions
library(xtable)

NREP <- 10000     # Monte Carlo replicates
n <- 30        # Sample size

set.seed(2022)  # Seed
# Generate model regressors
x0 <- rep(1,n)
x1 <- runif(n)
z0 <- rep(1,n)
z1 <- runif(n)

X <- cbind(x0,x1)  # Mean regressor matrix
Z <- cbind(z0,z1)  # Precision regressor matrix

# Matrices to store the criteria of the simulation
crit <- array(0, dim=c(NREP,10,16))

# Matrix to store the true values of the parameters
par <- matrix(0, 16, 4)

# Scenarios
par[1,] <- c(5, -2.63, 2.3, 2.8)
par[2,] <- c(5, -2.63, 5, 4)
par[3,] <- c(5, 7.24, 2.3, 2.8)
par[4,] <- c(5, 7.24, 5, 4)
par[5,] <- c(5, -2.63, 1.71, 1.55)
par[6,] <- c(5, -2.63, 3.23, 1.16)
par[7,] <- c(5, 7.24, 1.71, 1.55)
par[8,] <- c(5, 7.24, 3.23, 1.16)
par[9,] <- c(3.22, -1.5, 2.3, 2.8)
par[10,] <- c(3.22, -1.5, 5, 4)
par[11,] <- c(3.22, 1.79, 2.3, 2.8)
par[12,] <- c(3.22, 1.79, 5, 4)
par[13,] <- c(3.22, -1.5, 1.71, 1.55)
par[14,] <- c(3.22, -1.5, 3.23, 1.16)
par[15,] <- c(3.22, 1.79, 1.71, 1.55)
par[16,] <- c(3.22, 1.79, 3.23, 1.16)

for (k in 1:16) {
  cont <- 0 # convergence failure counts

  # Linear predictor vector
  eta1 <- X%*%par[k,1:2]
  eta2 <- Z%*%par[k,3:4]

  if(k >= 1 && k <= 4){
    mu <- as.vector((eta1)^2)
    phi <- as.vector((eta2)^2)
  }

  if(k >= 5 && k <= 8){
    mu <- as.vector((eta1)^2)
    phi <- as.vector(exp(eta2))
  }

  if(k >= 9 && k <= 12){
    mu <- as.vector(exp(eta1))
    phi <- as.vector((eta2)^2)
  }

  if(k >= 13 && k <= 16){
    mu <- as.vector(exp(eta1))
    phi <- as.vector(exp(eta2))
  }

  # Variance of response variable
  Vmu <-  mu*(1+mu)   # Variance function
  vary <-  Vmu/phi

  # Skewness of response variable
  ske <- (2*(1+phi)*(1+2*mu) / (phi-1)) * sqrt(phi/(Vmu*((1+phi)^2)))

  # Kurtosis of response variable
  kurt <- 6*((5*phi -1)/((phi-2)*(phi-1)) + (phi/(Vmu*(phi-2)*(phi-1))))

  cat("Scenario:", k, "\n")
  cat("Summary mu\n"); print(summary(mu))
  cat("Summary phi\n"); print(summary(phi))
  cat("Summary Variance of y\n"); print(summary(vary))
  cat("Summary Skewness of y\n"); print(summary(ske))
  cat("Summary kurtosis of y\n"); print(summary(kurt))

  # Monte Carlo Simulation
  i <- 1
  while(i <= NREP)
  {
    # Print replica ratio
    prop = (i/NREP)*100
    if(prop==0 || prop==10 || prop==20 || prop==30 || prop==40 || prop==50 || prop==60 || prop==70
       || prop==80|| prop== 90 || prop==100) cat(paste(prop,"%"),"\n")

    # Generate the response variable
    y <- rBP(n,mu,phi)

    # tryCatch to catch errors
    tryCatch({
    if(k >= 1 && k <= 4){
      # Estimate the model
      fit <- gamlss(y~x1,sigma.formula=~z1,family = BP(mu.link = "sqrt",sigma.link = "sqrt"), trace=FALSE)
    }

    if(k >= 5 && k <= 8){
      # Estimate the model
      fit <- gamlss(y~x1,sigma.formula=~z1,family = BP(mu.link = "sqrt",sigma.link = "log"), trace=FALSE)
    }

    if(k >= 9 && k <= 12){
      # Estimate the model
      fit <- gamlss(y~x1,sigma.formula=~z1,family = BP(mu.link = "log",sigma.link = "sqrt"), trace=FALSE)
    }

    if(k >= 13 && k <= 16){
      # Estimate the model
      fit <- gamlss(y~x1,sigma.formula=~z1,family = BP(mu.link = "log",sigma.link = "log"), trace=FALSE)
    }

      # if converge
      if(fit$converged == TRUE)
      {
        crit[i,,k] <- criteria(fit)
        i <- i + 1
      } else { # if it does not converge
        cont <- cont + 1
        print(c("Non-convergence",i,cont))
      }

    }, error = function(e) {
      # If an error occurs, print the message and repeat the loop
      print(paste("Error in the replica", i, ":", e$message))
      cont <- cont + 1
    })
  }
}

mean_crit_sqrt <- colMeans(crit[,,1:8])
result <- cbind(mean_crit_sqrt)
xtable(result, digits=3)

mean_crit_log <- colMeans(crit[,,9:16])
result <- cbind(mean_crit_log)
xtable(result, digits=3)
