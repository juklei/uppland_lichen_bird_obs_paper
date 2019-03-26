## bird bpo model with a random intercept
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
        nseen[i,k,y] ~ dbin(occ_true[i,k,y]*p_det[k,y], nvisits[i,k,y])
        nseen_sim[i,k,y] ~ dbin(occ_true[i,k,y]*p_det[k,y], nvisits[i,k,y])
      }
    logit(p_det[k,y]) <- alpha_p_det[k] + beta_obs_time*obs_time[y]
    }
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        occ_true[i,k,y] ~ dbern(p_occ[i,k,y])
        logit(p_occ[i,k,y]) <- alpha[k] + year_effect[k,y] + 
                               beta_ud[k]*understory_density[i,1] +
                               beta_cd[k]*canopy_density[i,1] +
                               beta_stand_dbh[k]*stand_dbh[i,1]
      }
    }
  }
  
  for(k in 1:nspecies){
    for(y in 1:nyears){
      year_effect[k,y] ~ dnorm(0, tau_year[k])
    }
  }
  
  ## Priors:
  
  for(k in 1:nspecies){
    alpha_p_det[k] ~ dnorm(0, 0.001)
    alpha[k] ~ dnorm(0, 0.001)
    beta_ud[k] ~ dnorm(0, 0.001)
    beta_cd[k] ~ dnorm(0, 0.001)
    beta_stand_dbh[k] ~ dnorm(0, 0.001)
    sigma_year[k] ~ dgamma(0.001, 0.001)
    tau_year[k] <- 1/sigma_year[k]^2
  }
  
  beta_obs_time ~ dnorm(0.5, 0.01)

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
        sq[i,k,y] <- (nseen[i,k,y] -
                      occ_true[i,k,y]*p_det[k,y]*nvisits[i,k,y])^2
        sq_sim[i,k,y] <- (nseen_sim[i,k,y] -
                          occ_true[i,k,y]*p_det[k,y]*nvisits[i,k,y])^2
      }
    }
  }

  fit <- sum(sq[,,])
  fit_sim <- sum(sq_sim[,,])
  p_fit <- step(fit_sim - fit)
  
  ## Predictions:
  
  # Understory density:
  for(m in 1:length(ud_pred)){
    for(k in 1:nspecies){
      occ_true_ud[m,k] ~ dbern(p_occ_ud[m,k])
      logit(p_occ_ud[m,k]) <- alpha[k] + beta_ud[k]*ud_pred[m]
    }
    r_ud[m] <- sum(p_occ_ud[m,])
  }

  ## Canopy density:
  for(m in 1:length(cd_pred)){
    for(k in 1:nspecies){
      occ_true_cd[m,k] ~ dbern(p_occ_cd[m,k])
      logit(p_occ_cd[m,k]) <- alpha[k] + beta_cd[k]*cd_pred[m]
    }
    r_cd[m] <- sum(p_occ_cd[m,])
  }
  
  ## Stand dbh:
  for(m in 1:length(stand_dbh_pred)){
    for(k in 1:nspecies){
      occ_true_stand_dbh[m,k] ~ dbern(p_occ_stand_dbh[m,k])
      logit(p_occ_stand_dbh[m,k]) <- alpha[k] + 
                                     beta_stand_dbh[k]*stand_dbh_pred[m]
    }
    r_stand_dbh[m] <- sum(p_occ_stand_dbh[m,])
  }
  
  
  ## Plot level richness:
  for(i in 1:nsites){
    for(y in 1:nyears){
      r_year[i,y] <- sum(occ_true[i,,y])
    }
    r_plot[i] <- sum(r_year[i,])/nyears
  }
  
}


