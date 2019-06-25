## bird and lichen richness modeled simultaniously with covariates
##
## First edit: 20190624
## Last edit: 20190624
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Ecological process model:
  for(p in 1:nsites){
    for(o in 1:2) {
    r_mean[p,o] ~ dnorm(richness[p,o], tau_richness[p,o])
    ## Moment matching:
    tau_richness[p,o] <- 1/r_sd[p,o]^2
    richness[p,o] <- alpha[o] + 
                     beta_ud[o]*ud[p,1] +
                     beta_cd[o]*cd[p,1] +
                     beta_sdbh[o]*stand_dbh[p,1]
  }}
  
  ## Priors:
  
  for(o in 1:2){
    alpha[o] ~ dnorm(0, 0.001)
    beta_ud[o] ~ dnorm(0, 0.001)
    beta_cd[o] ~ dnorm(0, 0.001)
    beta_sdbh[o] ~ dnorm(0, 0.001)
  }

  
  ## Predictions:
  
  # # Understory density:
  # for(m in 1:length(ud_pred)){
  #   r_ud[m] <- alpha + beta_ud*ud_pred[m]
  # }
  
}


