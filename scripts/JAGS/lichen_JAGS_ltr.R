## lichen ltr model
##
## First edit: 20190129
## Last edit: 20190131
##
## Author: Julian Klein, Matt Low

model{
  
  ## Likelihood:
  
  ## Tree level:
  for(i in 1:nobs){
  
    richness[i] ~ dbin(p, richness_true[i])
    
    richness_true[i] ~ dpois(lambda[i])
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
  p ~ dnorm(0.8, 16)
  beta_pine ~ dnorm(0, 0.01)
  beta_spruce ~ dnorm(0, 0.01)
  beta_stem_dbh ~ dnorm(0.25, 0.01)
  
  sigma_plot ~ dunif(0, 5)
  tau_plot <- 1/sigma_plot^2
  
  alpha_plot_mean ~ dnorm(2, 0.04)
  beta_stand_dbh ~ dnorm(0.25, 0.01)
  beta_cdens ~ dnorm(0, 0.01)
  beta_udens ~ dnorm(0, 0.01)
  
  # sigma_block ~ dunif(0, 5)
  # tau_block <- 1/sigma_block^2

  ## Predictions:
  for(m in 1:length(ud_pred)){
    
    log(out_d[m]) <- alpha_plot_mean + beta_udens*ud_pred[m]
    log(out_p[m]) <- alpha_plot_mean + beta_udens*ud_pred[m] + beta_pine
    log(out_s[m]) <- alpha_plot_mean + beta_udens*ud_pred[m] + beta_spruce
    mean_out[m] <- (out_d[m] + out_p[m] + out_s[m])/3
    
  }
   
}

