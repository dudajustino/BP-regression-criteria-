
criteria <- function(model){
  y <- model$y
  n <- model$N
  p <- length(model$mu.coefficients)
  q <- length(model$sigma.coefficients)
  r <- p + q
  mu <-  model$mu.fv
  sigma <-  model$sigma.fv
  eta <- model$mu.lp
  mu.link <- model$mu.link
  sigma.link <- model$sigma.link
  influence <- influence.BP(model)
  hat <- influence$hat
  ychap <- influence$ychap

  ############ PRESS #################

  press <- function(res) sum((res/(1-hat))^2)

  press_quantile <- press(residuals.BP(model, type="quantile"))
  press_sweighted1 <- press(residuals.BP(model, type="sweighted1"))
  press_sweighted2 <- press(residuals.BP(model, type="sweighted2"))

  ############ P2 #################

  sst <- sum((ychap-mean(ychap))^2)

  p2_quantile <- 1 - press_quantile/((n/(n-r))^2*sst)
  p2_sweighted1 <- 1 - press_sweighted1/((n/(n-r))^2*sst)
  p2_sweighted2 <- 1 - press_sweighted2/((n/(n-r))^2*sst)

  ###### P2 - correction (Espinheira et al., 2019) ######

  p2_quantile_c <- 1-(1-p2_quantile)*((n-1)/(n-r))
  p2_sweighted1_c <- 1-(1-p2_sweighted1)*((n-1)/(n-r))
  p2_sweighted2_c <- 1-(1-p2_sweighted2)*((n-1)/(n-r))

  ######### R2 likelihood ratio #################

  model0 <- gamlss(y~1, sigma.formula=~1, family = BP(mu.link = mu.link, sigma.link = sigma.link),
                   control = model$control, trace=FALSE)
  mu0 <-  model0$mu.fv
  sigma0 <-  model0$sigma.fv
  L0 <-  exp(log.like(y,mu0,sigma0))
  L1 <-  exp(log.like(y,mu,sigma))

  r2_RV <- 1-(L0/L1)^(2/n)

  ######### R2 likelihood ratio - correction (Bayer and Cribari-Neto, 2017) #################

  alpha <- 0.4
  delta <- 1

  r2_RV_c <- 1-(1-r2_RV)*((n-1)/(n-(1+alpha)*p-(1-alpha)*q))^delta

  ######### R2 Ferrari and Cribari-Neto (2004) #################

  if(mu.link == "log"){
    g <- log(y)
  } else if(mu.link=="sqrt"){
    g <- sqrt(y)
  } else {
    g <- y
  }

  r2_FC <- cor(eta,g)^2

  ###### R2 Ferrari and Cribari-Neto - correction (Bayer and Cribari-Neto, 2017) ######

  r2_FC_c <- 1-(1-r2_FC)*((n-1)/(n-r))

  medidas = c(p2_quantile, p2_quantile_c, p2_sweighted1, p2_sweighted1_c, p2_sweighted2, p2_sweighted2_c,
              r2_FC, r2_FC_c, r2_RV, r2_RV_c)

  return(medidas)
}
