## bird bpo model
##
## First edit: 20190307
## Last edit: 20190311
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observational model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        nseen[i,k,y] ~ dbin(occ_true[i,k,y]*p_det[k], nvisits[i,k,y])
      }
    }
  }
  
  ## Ecological process model:
    for(k in 1:nspecies){
      for(y in 1:nyears){
        for(i in 1:nsites){
          occ_true[i,k,y] ~ dbern(p_occ[i,k,y])
          logit(p_occ[i,k,y]) <- alpha[k] + beta_ud[k]*understory_density[i,1]
      }
    }
  }
  
  # ## Year effect:
  # for(k in 1:nspecies){
  #   for (y in 1:nyears){
  #     alpha_year_mean[k,y] ~ dnorm(mu_year, tau_year)
  #   }
  # }
  
  ## Priors:
  
  for(k in 1:nspecies){
    p_det[k] ~ dunif(0, 1)
    alpha[k] ~ dnorm(0, 0.01)
    beta_ud[k] ~ dnorm(0, 0.01)
  }

  ## Model validation:
  
  #...
  
  ## Predictions:
  
  for(m in 1:length(ud_pred)) {
    for(k in 1:nspecies){
      occ_true_pred[m,k] ~ dbern(p_occ_pred[m,k])
      logit(p_occ_pred[m,k]) <- alpha[k] + beta_ud[k]*ud_pred[m]
    }
    richness[m] <- sum(occ_true_pred[m,])
  }
  
}


