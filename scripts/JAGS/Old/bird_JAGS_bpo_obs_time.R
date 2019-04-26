## bird bpo model with a fixed intercept
##
## First edit: 20190307
## Last edit: 20190314
##
## Author: Julian Klein, Matt Low

model{
  
  ## Likelihood:
  
  ## Observational model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        obs[i,k,y] ~ dbern(occ_true[i,k,y]*p_det[i,k,y])
        p_det[i,k,y] <- obs_time[i,k,y]/(param_obs_time + obs_time[i,k,y])
      }
    }
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        occ_true[i,k,y] ~ dbern(p_occ)
      }
    }
  }
  
  ## Priors:
  
  param_obs_time ~ dunif(0.05, 3)
  p_occ ~ dunif(0,1)

  ## Plot level richness:
  for(i in 1:nsites){
    for(y in 1:nyears){
      r_year[i,y] <- sum(occ_true[i,,y])
    }
  }
  
}


