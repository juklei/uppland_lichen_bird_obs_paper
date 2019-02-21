## lichen ltr model
##
## First edit: 20190129
## Last edit: 20190220
##
## Author: Julian Klein, Matt Low

model{
  
  ## Likelihood:
  
  ## Tree level:
  for(i in 1:nobs){
  
    richness[i] ~ dbin(p[i], richness_true[i])
    richness_sim[i] ~ dbin(p[i], richness_true[i])
    
    p[i] ~ dbeta(a[i], b[i])
    ## Moment matching:
    a[i] <- (mu_p^2 - mu_p^3 - mu_p*sigma_p^2)/sigma_p^2
    b[i] <- (mu_p - 2*mu_p^2 + mu_p^3 - sigma_p^2 + mu_p*sigma_p^2)/sigma_p^2
    
    richness_true[i] ~ dpois(lambda[i]) #T(1, 50)
    log(lambda[i]) <- alpha_plot[plot[i]] + 
                      beta_pine*tree_sp_pine[i] + 
                      beta_spruce*tree_sp_spruce[i] + 
                      beta_stem_dbh*stem_dbh[i, 1]
    
  }
  
  ## Plot level:
  for(j in 1:nplot){

    alpha_plot[j] ~ dnorm(mu[j], tau_plot)
    mu[j] <- alpha_plot_mean +
             beta_stand_dbh*stand_dbh[j, 1] +
             beta_cdens*canopy_density[j, 1] +
             beta_udens*understory_density[j, 1]

  }

  ## Block level:
  # for (k in 1:nblock){
  # 
  #   alpha_plot_mean[k] ~ dnorm(mu2, tau_block)
  # 
  # }
  
  ## Priors:
  
  # p ~ dbeta(7.2, .8)
  # alpha_p ~ dnorm(-0.05, 0.01)
  # beta_p ~ dnorm(-0.05, 0.01)
  sigma_p ~ dunif(0, 1)
  
  beta_pine ~ dnorm(0, 0.01)
  beta_spruce ~ dnorm(0, 0.01)
  beta_stem_dbh ~ dnorm(0.25, 0.01)
  
  alpha_plot_mean ~ dnorm(2, 0.04)
  beta_stand_dbh ~ dnorm(0.25, 0.01)
  beta_cdens ~ dnorm(0, 0.01)
  beta_udens ~ dnorm(0, 0.01)
  sigma_plot ~ dunif(0, 5)
  tau_plot <- 1/sigma_plot^2
  
  # sigma_block ~ dunif(0, 5)
  # tau_block <- 1/sigma_block^2

  ## Model validation:

  ## Bayesian p-value:
  mean_richness <- mean(richness[])
  mean_richness_sim <- mean(richness_sim[])
  p_mean <- step(mean_richness_sim - mean_richness)

  ## Coefficient of variation:
  cv_richness <- sd(richness[])/mean_richness
  cv_richness_sim <- sd(richness_sim[])/mean_richness_sim
  p_cv <- step(cv_richness - cv_richness_sim)

  ## Model fit:
  for(m in 1:nobs){

    sq[m] <- (richness[m] - p[m]*richness_true[m])^2
    sq_sim[m] <- (richness_sim[m] - p[m]*richness_true[m])^2
    
  }

  fit <- sum(sq[])
  fit_sim <- sum(sq_sim[])
  p_fit <- step(fit_sim - fit)

  #R2 <- 1 - (fit/sum(sq_mean[]))
  #RMSE <- sqrt(mean(sq[]^2))
  
  ## Predictions:
  for(n in 1:length(ud_pred)){

    log(out_d[n]) <- alpha_plot_mean + beta_udens*ud_pred[n]
    log(out_p[n]) <- alpha_plot_mean + beta_udens*ud_pred[n] + beta_pine
    log(out_s[n]) <- alpha_plot_mean + beta_udens*ud_pred[n] + beta_spruce
    mean_out[n] <- (out_d[n] + out_p[n] + out_s[n])/3

  }
   
}

