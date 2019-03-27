## bird and lichen richness modeled simultaniously with covariates
##
## First edit: 20190326
## Last edit: 20190326
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Ecological process model:
  for(i in 1:nobs){
    r_mean[i] ~ dnorm(richness[i], tau_richness[i])
    ## Moment matching:
    tau_richness[i] <- 1/r_sd[i]^2
    richness[i] <- alpha + 
                   site_effect[sites[i]] + 
                   year_effect[obs_year[i]] + 
                   beta_ud*ud[i,1]
  }
  
  for(j in 1:nsites){
    site_effect[j] ~ dnorm(0, tau_site)
  }
  
  for(y in years){
    year_effect[y] ~ dnorm(0, tau_year)
  }
  
  ## Priors:
  
  alpha ~ dnorm(10, 0.001)
  beta_ud ~ dnorm(0.5, 0.001)
  
  sigma_site ~ dgamma(0.001, 0.001)
  tau_site <- 1/sigma_site^2
  
  sigma_year ~ dgamma(0.001, 0.001)
  tau_year <- 1/sigma_year^2

  ## Predictions:
  
  # Understory density:
  for(m in 1:length(ud_pred)){
    r_ud[m] <- alpha + beta_ud*ud_pred[m]
  }
  
}


