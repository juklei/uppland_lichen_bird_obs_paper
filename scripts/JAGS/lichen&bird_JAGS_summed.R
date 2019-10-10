## bird and lichen richness modeled simultaniously with covariates
##
## First edit: 20190716
## Last edit: 20191007
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Ecological process model:
  for(i in 1:nobs){
    r_mean[i] ~ dnorm(richness[i], 1/r_sd[i]^2)
    # r_sim[i] ~ dnorm(richness[i], 1/r_sd[i]^2)
    richness[i] <- alpha + residual[i] +
                   beta_ud*ud[i,1] +
                   beta_od*od[i,1] +
                   beta_sdbh*sdbh[i,1]
    residual[i] ~ dnorm(0, 1/plot_sd^2)
  }
  
  ## Priors:
  
  alpha ~ dnorm(0, 0.001)
  beta_ud ~ dnorm(0, 0.001)
  beta_od ~ dnorm(0, 0.001)
  beta_sdbh ~ dnorm(0, 0.001)
  plot_sd ~ dunif(0.0001, 10)

#   ## Model validation:
#   
#   ## Bayesian p-value:
#   mean_r_mean <- mean(r_mean[])
#   mean_r_sim <- mean(r_sim[])
#   p_mean <- step(mean_r_sim - mean_r_mean)
# 
#   ## Coefficient of variation:
#   cv_r_mean <- sd(r_mean[])/mean_r_mean
#   cv_r_sim <- sd(r_sim[])/mean_r_sim
#   p_cv <- step(cv_r_sim - cv_r_mean)
# 
#   ## Model fit:
#   for(m in 1:nobs){
#     sq[m] <- (r_mean[m] - richness[m])^2
#     sq_sim[m] <- (r_sim[m] - richness[m])^2
#   }
# 
#   fit <- sum(sq[])
#   fit_sim <- sum(sq_sim[])
#   p_fit <- step(fit_sim - fit)
#   
# }

  ## Predictions:
  
  # Understory density:
  for(m in 1:length(ud_pred)){
    r_ud[m] <- alpha + beta_ud*ud_pred[m]
  }

  # Overstory density:
  for(m in 1:length(od_pred)){
    r_od[m] <- alpha + beta_od*od_pred[m]
  }
  
  # # Sdbh:
  # for(m in 1:length(height_pred)){
  #   r_height[m] <- alpha + beta_height*height_pred[m]
  # }
  
}


