## bird bpr model
##
## First edit: 20190328
## Last edit: 20190617
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
    log(lambda[i])  <- plot_richness[plot[i]] + year_effect[obs_year[i]]
  }
  
  ## Year effect:
  for(y in years){
    year_effect[y] ~ dnorm(0, tau_year)
  }
  
  ## Priors:
  
  param_obs_time ~ dunif(1, 60)
  
  for(p in 1:nsites){
    plot_richness[p] ~ dgamma(0.001, 0.001)
  }
  
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
  
}


