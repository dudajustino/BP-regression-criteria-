
# Clear all existing objects from the workspace
rm(list=ls())

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R")
source("criteria.R")

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions
library(xtable)

NREP <- 10000     # Monte Carlo replicates
n <- 40        # Sample size

set.seed(2022)  # Seed
# Generate model regressors
x0 <- rep(1,n)
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- runif(n)
z0 <- rep(1,n)
z1 <- runif(n)
z2 <- runif(n)

X <- cbind(x0,x1,x2,x3,x4) # mean regressor matrix
Z <- cbind(z0,z1,z2) # Precision regressor matrix

# Matrices to store the criteria of the simulation
crit <- array(0, dim=c(NREP,12,8))
select_correct <- matrix(0, nrow=12, ncol=2)

# Matrix to store the true values of the parameters
par <- matrix(0, 2, 8)

par[1,] <- c(1.1, 0.8, 1.3, 1.2, 1.3, 3, -1.0, 0.9)
par[2,] <- c(1.1, 0.8, 1.3, 1.2, 1.3, 4.5, -1.2, 0.7)

for (k in 1:2) {

  cont <- 0 # convergence failure counts

  # Linear predictor vector
  eta1 <- X%*%par[k,1:5]
  eta2 <- Z%*%par[k,6:8]

  mu <- as.vector(exp(eta1))
  phi <- as.vector(exp(eta2))

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
  valid <- 0
  i <- 1
  while(i <= NREP)
  {
    # Print replica ratio
    prop = (i/NREP)*100
    if(prop==0 || prop==10 || prop==20 || prop==30 || prop==40 || prop==50 || prop==60 || prop==70
       || prop==80|| prop== 90 || prop==100) cat(paste(prop,"%"),"\n")

    # Generate the response variable
    y <- rBP(n,mu,phi)

    tryCatch({
      # correct model
      fit4 <- gamlss(y~x1+x2+x3+x4,sigma.formula=~z1+z2,family = BP(mu.link = "log",sigma.link = "log"), trace=FALSE)

      # covariate omission specificatio
      fit1 <-  gamlss(y~x1,sigma.formula=~z1+z2,family = BP(mu.link = "log",sigma.link = "log"), trace=FALSE)
      fit2 <-  gamlss(y~x1+x2,sigma.formula=~z1+z2,family = BP(mu.link = "log",sigma.link = "log"), trace=FALSE)
      fit3 <-  gamlss(y~x1+x2+x3,sigma.formula=~z1+z2,family = BP(mu.link = "log",sigma.link = "log"), trace=FALSE)

      # if converge
      if(fit1$converged == TRUE && fit2$converged == TRUE && fit3$converged == TRUE && fit4$converged == TRUE)
      {
        c1 <- c(criteria(fit1), AIC(fit1), BIC(fit1))
        c2 <- c(criteria(fit2), AIC(fit2), BIC(fit2))
        c3 <- c(criteria(fit3), AIC(fit3), BIC(fit3))
        c4 <- c(criteria(fit4), AIC(fit4), BIC(fit4))

        C <- cbind(c1,c2,c3,c4)

        pos <- (k-1)*4

        crit[i,,pos+1] <- c1
        crit[i,,pos+2] <- c2
        crit[i,,pos+3] <- c3
        crit[i,,pos+4] <- c4

        for(j in 1:12){
          if(j <= 10){
            escolhido <- which.max(C[j,])
          } else {
            escolhido <- which.min(C[j,])
          }

        if(escolhido == 4){
          select_correct[j,k] <- select_correct[j,k] + 1
        }
      }

      i <- i + 1
      valid <- valid + 1

      } else{ # if it does not converge
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

mean_criterios <- colMeans(crit[,,1:4])
result <- cbind(mean_criterios)
xtable(result, digits=3)

mean_criterios <- colMeans(crit[,,5:8])
result <- cbind(mean_criterios)
xtable(result, digits=3)

prop_correct <- select_correct / NREP * 100
rownames(prop_correct) <- c(
  "P2_q", "P2_q_c", "P2_sw1", "P2_sw1_c", "P2_sw2", "P2_sw2_c",
  "R2_FC", "R2_FC_c", "R2_RV", "R2_RV_c",
  "AIC", "BIC"
)
colnames(prop_correct) <- c("C1","C2")
xtable(prop_correct, digits=3)
