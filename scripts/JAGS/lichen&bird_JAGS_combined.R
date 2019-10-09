## bird and lichen richness modeled simultaniously with covariates
##
## First edit: 20190624
## Last edit: 201901004
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Ecological process model:
  for(i in 1:nobs){
    r_mean[i] ~ dnorm(richness[i], 1/r_sd[i]^2)
    # r_sim[i] ~ dnorm(richness[i], 1/r_sd[i]^2)
    richness[i] <- alpha + residual[i] +
                   beta_org*org[i] +
                   beta_ud*ud[i,1] +
                   beta_od*od[i,1] +
                   beta_sdbh*sdbh[i,1] +
                   int_ud*org[i]*ud[i,1] +
                   int_od*org[i]*od[i,1] +
                   int_sdbh*org[i]*sdbh[i,1]
    residual[i] ~ dnorm(0, 1/plot_sd^2)
  }
  
  ## Priors:
  
  alpha ~ dnorm(0, 0.001)
  beta_org ~ dnorm(0, 0.001)
  beta_ud ~ dnorm(0, 0.001)
  beta_od ~ dnorm(0, 0.001)
  beta_sdbh ~ dnorm(0, 0.001)
  int_ud ~ dnorm(0, 0.001)
  int_od ~ dnorm(0, 0.001)
  int_sdbh ~ dnorm(0, 0.001)
  plot_sd ~ dgamma(0.001, 0.001)

  # ## Model validation:
  # 
  # ## Bayesian p-value:
  # mean_r_mean <- mean(r_mean[])
  # mean_r_sim <- mean(r_sim[])
  # p_mean <- step(mean_r_sim - mean_r_mean)
  # 
  # ## Coefficient of variation:
  # cv_r_mean <- sd(r_mean[])/mean_r_mean
  # cv_r_sim <- sd(r_sim[])/mean_r_sim
  # p_cv <- step(cv_r_mean - cv_r_sim)
  # 
  # ## Model fit:
  # for(m in 1:nobs){
  #   sq[m] <- (r_mean[m] - richness[m])^2
  #   sq_sim[m] <- (r_sim[m] - richness[m])^2
  # }
  # 
  # fit <- sum(sq[])
  # fit_sim <- sum(sq_sim[])
  # p_fit <- step(fit_sim - fit)
  
  ## Predictions:
  
  # Understory density:
  for(m in 1:length(ud_pred)){
    rb_ud[m] <- alpha + beta_ud*ud_pred[m]
    rl_ud[m] <- alpha + beta_org + (beta_ud + int_ud)*ud_pred[m]
  }

  # Overstory density:
  for(m in 1:length(od_pred)){
    rb_od[m] <- alpha + beta_od*od_pred[m]
    rl_od[m] <- alpha + beta_org + (beta_od + int_od)*od_pred[m]
  }

  # # Sdbh:
  # for(m in 1:length(sdbh_pred)){
  #   rb_sdbh[m] <- alpha + beta_sdbh*sdbh_pred[m]
  #   rl_sdbh[m] <- alpha + beta_org + (beta_sdbh + int_sdbh)*sdbh_pred[m]
  # }
  
}


