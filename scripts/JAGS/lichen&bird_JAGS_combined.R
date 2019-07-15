## bird and lichen richness modeled simultaniously with covariates
##
## First edit: 20190624
## Last edit: 20190715
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
                   beta_org*org[i] +
                   beta_ud*ud[i,1] +
                   beta_cd*cd[i,1] +
                   beta_sdbh*stand_dbh[i,1] +
                   int_ud*org[i]*ud[i,1] +
                   int_cd*org[i]*cd[i,1] +
                   int_sdbh*org[i]*stand_dbh[i,1]
  }
  
  ## Priors:
  
  alpha ~ dnorm(0, 0.001)
  beta_org ~ dnorm(0, 0.001)
  beta_ud ~ dnorm(0, 0.001)
  beta_cd ~ dnorm(0, 0.001)
  beta_sdbh ~ dnorm(0, 0.001)
  int_ud ~ dnorm(0, 0.001)
  int_cd ~ dnorm(0, 0.001)
  int_sdbh ~ dnorm(0, 0.001)

  ## Predictions:
  
  # Understory density:
  for(m in 1:length(ud_pred)){
    rb_ud[m] <- beta_ud*ud_pred[m]
    rl_ud[m] <- (beta_ud + int_ud)*ud_pred[m]
  }

  # Canopyy density:
  for(m in 1:length(cd_pred)){
    rb_cd[m] <- beta_cd*cd_pred[m]
    rl_cd[m] <- (beta_cd + int_cd)*cd_pred[m]
  }

  # Stand age:
  for(m in 1:length(sdbh_pred)){
    rb_sdbh[m] <- alpha + beta_sdbh*sdbh_pred[m]
    rl_sdbh[m] <- alpha + beta_org + (beta_sdbh + int_sdbh)*sdbh_pred[m]
  }
  
}


