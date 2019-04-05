## bird bpr model
##
## First edit: 20190328
## Last edit: 20190328
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observational level:
  for(i in 1:nobs){
    richness[i] ~ dbin(p_det[i], richness_true[i])
    # richness_sim[i] ~ dbin(p_det[i], richness_true[i])
    p_det[i] <- obs_time[i]/(param_obs_time + obs_time[i])
  }
  
  ## Process model level:
  for(i in 1:nobs){
    richness_true[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha_plot_mean + 
                      plot_effect[plot[i]] + 
                      year_effect[obs_year[i]] +
                      beta_stand_dbh*stand_dbh[i,1] +
                      beta_cd*cd[i,1] + 
                      beta_ud*ud[i,1] #+
                      # beta_spruce*spruce[i,1] +
                      # beta_pine*pine[i,1] +
                      # beta_umbr*umbrella[i,1] +
                      # beta_dec*dec[i,1]
  }
  
  ## Group effects:
  
  for(j in 1:nsites){
    plot_effect[j] ~ dnorm(0, tau_plot)
  }
  
  for(y in years){
    year_effect[y] ~ dnorm(0, tau_year)
  }
  
  ## Priors:
  
  param_obs_time ~ dunif(0.05, 3)
  alpha_plot_mean ~ dnorm(2, 0.01)
  beta_stand_dbh ~ dnorm(0, 0.001)
  beta_cd ~ dnorm(0, 0.001)
  beta_ud ~ dnorm(0, 0.001)
  # beta_spruce ~ dnorm(0, 0.001)
  # beta_pine ~ dnorm(0, 0.001)
  # beta_dec ~ dnorm(0, 0.001)
  # beta_umbr ~ dnorm(0, 00.1)
  
  sigma_plot ~ dgamma(0.001, 0.001)
  tau_plot <- 1/sigma_plot^2

  sigma_year ~ dgamma(0.001, 0.001)
  tau_year <- 1/sigma_year^2
  
  # ## Model validation:
  # 
  # ## Bayesian p-value:
  # mean_richness <- mean(richness[])
  # mean_richness_sim <- mean(richness_sim[])
  # p_mean <- step(mean_richness_sim - mean_richness)
  # 
  # ## Coefficient of variation:
  # cv_richness <- sd(richness[])/mean_richness
  # cv_richness_sim <- sd(richness_sim[])/mean_richness_sim
  # p_cv <- step(cv_richness - cv_richness_sim)
  # 
  # ## Model fit:
  # for(m in 1:nobs){
  #   sq[m] <- (richness[m] - p_det[m]*richness_true[m])^2
  #   sq_sim[m] <- (richness_sim[m] - p_det[m]*richness_true[m])^2
  # }
  # 
  # fit <- sum(sq[])
  # fit_sim <- sum(sq_sim[])
  # p_fit <- step(fit_sim - fit)
  
  ## Predictions:
  
  # Understory density:
  for(m in 1:length(ud_pred)){
    log(r_ud[m]) <- alpha_plot_mean + beta_ud*ud_pred[m]
  }
  
  # Canopy density:
  for(m in 1:length(cd_pred)){
    log(r_cd[m]) <- alpha_plot_mean + beta_cd*cd_pred[m]
  }
  
  # Stand_dbh density:
  for(m in 1:length(stand_dbh_pred)){
    log(r_stand_dbh[m]) <- alpha_plot_mean + beta_stand_dbh*stand_dbh_pred[m]
  }
  
}


