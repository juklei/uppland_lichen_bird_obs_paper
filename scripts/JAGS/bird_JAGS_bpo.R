## bird bpo model
##
## First edit: 20190307
## Last edit: 20190307
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observational model:
  for(k in 1:nspecies){
    for(i in 1:nsites){
      nseen[i,k] ~ dbin(occ_true[i,k]*p_det[k], nvisits[i,k])
    }
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(i in 1:nsites){
      occ_true[i,k] ~ dbern(p_occ[k])
      #logit(p_occ[k]) <- alpha + beta_ud*understorey_density[i]
    }
  }
  
  ## Block level:
  # for (k in 1:nblock){
  #   alpha_plot_mean[k] ~ dnorm(mu2, tau_block)
  # }
  
  ## Priors:
  
  for(k in 1:nspecies){
    p_det[k] ~ dunif(0, 1)
    p_occ[k] ~ dunif(0, 1)
  }

  ## Model validation:
  
  #...
  
  ## Predictions:

  for(i in 1:nsites){
    richness[i] <- sum(occ_true[i,])
  }
  
}


