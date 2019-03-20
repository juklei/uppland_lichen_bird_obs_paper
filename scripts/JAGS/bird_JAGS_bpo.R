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
    logit(p_det[k,y]) <- alpha_p_det[k] + beta_obs_time*log(obs_time[k,y])
    }
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    for(y in 1:nyears){
      for(i in 1:nsites){
        occ_true[i,k,y] ~ dbern(p_occ[i,k,y])
        logit(p_occ[i,k,y]) <- alpha[k,y] + 
                               beta_ud[k]*understory_density[i,1] +
                               beta_cd[k]*canopy_density[i,1] +
                               beta_stand_dbh[k]*stand_dbh[i,1]
      }
    }
  }
  
  ## Species effect:
  for(k in 1:nspecies){
    for(y in 1:nyears){
     alpha[k,y] ~ dnorm(alpha_mean, tau_alpha)
    }
  }
  
  ## Priors:
  
  for(k in 1:nspecies){
    alpha_p_det[k] ~ dnorm(0, 0.0001)
    beta_ud[k] ~ dnorm(0, 0.0001)
    beta_cd[k] ~ dnorm(0, 0.0001)
    beta_stand_dbh[k] ~ dnorm(0, 0.0001)
  }
  
  beta_obs_time ~ dnorm(0.5, 0.001)
  alpha_mean ~ dnorm(0, 0.0001)
  sigma_alpha ~ dgamma(0.0001, 0.0001)
  tau_alpha <- 1/sigma_alpha^2

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
  
  ## Understory density:
  for(m in 1:length(ud_pred)){
    for(k in 1:nspecies){
      for(y in 1:nyears){
        occ_true_ud[m,k,y] ~ dbern(p_occ_ud[m,k,y])
        logit(p_occ_ud[m,k,y]) <- alpha[k,y] + beta_ud[k]*ud_pred[m]
      }
      ots_ud[m,k] <- sum(occ_true_ud[m,k,])/nyears
    }
    r_ud[m] <- sum(ots_ud[m,])
  }
  
  ## Canopy density:
  for(m in 1:length(cd_pred)){
    for(k in 1:nspecies){
      for(y in 1:nyears){
        occ_true_cd[m,k,y] ~ dbern(p_occ_cd[m,k,y])
        logit(p_occ_cd[m,k,y]) <- alpha[k,y] + beta_cd[k]*cd_pred[m]
      }
      ots_cd[m,k] <- sum(occ_true_cd[m,k,])/nyears
    }
    r_cd[m] <- sum(ots_cd[m,])
  }
  
  ## Stand dbh:
  for(m in 1:length(stand_dbh_pred)){
    for(k in 1:nspecies){
      for(y in 1:nyears){
        occ_true_stand_dbh[m,k,y] ~ dbern(p_occ_stand_dbh[m,k,y])
        logit(p_occ_stand_dbh[m,k,y]) <- alpha[k,y] + 
                                         beta_stand_dbh[k]*stand_dbh_pred[m]
      }
      ots_stand_dbh[m,k] <- sum(occ_true_stand_dbh[m,k,])/nyears
    }
    r_stand_dbh[m] <- sum(ots_stand_dbh[m,])
  }
  
  ## Plot level richness:
  for(i in 1:nsites){
    for(y in 1:nyears){
      r_year[i,y] <- sum(occ_true[i,,y])
    }
    r_plot[i] <- sum(r_year[i,])/nyears
  }
  
}


