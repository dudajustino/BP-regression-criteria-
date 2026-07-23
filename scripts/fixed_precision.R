
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
n <- 60        # Sample size

set.seed(2022)  # Seed
# Generate model regressors
x0 <- rep(1,n)
x1 <- runif(n)

PHI <- c(20, 50, 150)

X <- cbind(x0,x1)  # Mean regressor matrix

# Matrices to store the criteria of the simulation
crit <- array(0, dim=c(NREP,10,12))

# Matrix to store the true values of the parameters
par <- matrix(0, 4, 2)

# Four scenarios
par[1,] <- c(5, -2.63)    #g(\mu) = \sqrt
par[2,] <- c(5, 7.24)     #g(\mu) = \sqrt
par[3,] <- c(3.22, -1.5)  #g(\mu) = \log
par[4,] <- c(3.22, 1.79)  #g(\mu) = \log

for(l in 1:3){
  cont <- 0 # convergence failure counts

  phi1 <- PHI[l] #precision parameter

  for (k in 1:4) {
    # Linear predictor vector
    eta1 <- X%*%par[k,1:2]

    if(k >= 1 && k <= 2){
      mu <- as.vector((eta1)^2)
      phi <- as.vector(phi1)
    }

    if(k >= 3 && k <= 4){
      mu <- as.vector(exp(eta1))
      phi <- as.vector(phi1)
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

      if(k >= 1 && k <= 2){
        # Estimate the model
        fit <-  gamlss(y~x1,family = BP(mu.link = "sqrt", sigma.link = "identity"), trace=FALSE)
      }

      if(k >= 3 && k <= 4){
        # Estimate the model
        fit <-  gamlss(y~x1,family = BP(mu.link = "log", sigma.link = "identity"), trace=FALSE)
      }

      # if converge
      if(fit$converged == TRUE)
      {
        crit[i,,(k*3)-(3-l)] <- criteria(fit)

        i <- i + 1

      }else{ # if it does not converge
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
}

mean_crit_sqrt <- colMeans(crit[,,1:6])
result <- cbind(mean_crit_sqrt)
xtable(result, digits=3)

mean_crit_log <- colMeans(crit[,,7:12])
result <- cbind(mean_crit_log)
xtable(result, digits=3)
