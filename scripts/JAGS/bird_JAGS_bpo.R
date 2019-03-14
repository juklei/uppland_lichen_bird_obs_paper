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
        nseen_sim[i,k,y] ~ dbin(occ_true[i,k,y]*p_det[k], nvisits[i,k,y])
      }
    }
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        occ_true[i,k,y] ~ dbern(p_occ[i,k,y])
        logit(p_occ[i,k,y]) <- alpha_mean + 
                               spec_effect[k] + 
                               year_effect[y] + 
                               beta_ud[k]*understory_density[i,1]
      }
    }
  }
  
  ## Species effect:
  for(k in 1:nspecies){
    spec_effect[k] ~ dnorm(0, tau_spec)
  }
  
  ## Year effect:
  for (y in 1:nyears){
    year_effect[y] ~ dnorm(0, tau_year)
  }
  
  ## Priors:
  
  for(k in 1:nspecies){
    p_det[k] ~ dunif(0, 1)
    beta_ud[k] ~ dnorm(0, 0.01)
  }
  
  alpha_mean ~ dnorm(0, 0.001)
  sigma_spec ~ dgamma(0.001, 0.001)
  tau_spec <- 1/sigma_spec^2
  sigma_year ~ dgamma(0.001, 0.001)
  tau_year <- 1/sigma_year^2

  ## Model validation:
  
  ## Bayesian p-value:
  mean_nseen <- mean(nseen[,,])
  mean_nseen_sim <- mean(nseen_sim[,,])
  p_mean <- step(mean_nseen_sim - mean_nseen)
  
  ## Coefficient of variation:
  cv_nseen <- sd(nseen[,,])/mean_nseen
  cv_nseen_sim <- sd(nseen_sim[,,])/mean_nseen_sim
  p_cv <- step(cv_nseen - cv_nseen_sim)
  
  ## Model fit:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        sq[i,k,y] <- (nseen[i,k,y] - occ_true[i,k,y]*p_det[k]*nvisits[i,k,y])^2
        sq_sim[i,k,y] <- (nseen_sim[i,k,y] - 
                          occ_true[i,k,y]*p_det[k]*nvisits[i,k,y])^2
      }
    }
  }
  
  fit <- sum(sq[,,])
  fit_sim <- sum(sq_sim[,,])
  p_fit <- step(fit_sim - fit)
  
  ## Predictions:
  
  ## Understory density:
  for(m in 1:length(ud_pred)){
    for(k in 1:nspecies){
      occ_true_pred[m,k] ~ dbern(p_occ_pred[m,k])
      logit(p_occ_pred[m,k]) <- alpha_mean + beta_ud[k]*ud_pred[m]
    }
    richness[m] <- sum(occ_true_pred[m,])
  }
  
}


