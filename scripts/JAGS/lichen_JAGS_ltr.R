## lichen ltr model
##
## First edit: 20190129
## Last edit: 20190129
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Tree level:
  for(i in 1:nobs){
  
    richness[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha_plot[plot[i]] + 
                      beta_pine*tree_sp_pine[i] + 
                      beta_spruce*tree_sp_spruce[i] + 
                      beta_stem_dbh*stem_dbh[i, 1]
    
  }
  
  ## Plot level:
  for(j in 1:nplot){
  
    alpha_plot[j] ~ dgamma(mu[j]^2/sigma_plot^2, mu[j]/sigma_plot^2)
    log(mu[j]) <- alpha_plot_mean + block_effect[block[j]] +
                  beta_stand_dbh*stand_dbh[j, 1] +
                  beta_cdens*canopy_density[j, 1] +
                  beta_udens*understory_density[j, 1]
  
  } 

  ## Block level:
  for (k in 1:nblock){
  
    block_effect[k] ~ dnorm(0, tau_block)
  
  }
  
  ## Priors:
  
  beta_pine ~ dnorm(0, 0.001)
  beta_spruce ~ dnorm(0, 0.001)
  beta_stem_dbh ~ dnorm(0, 0.001)
  
  sigma_plot ~ dgamma(0.001, 0.001)
  alpha_plot_mean ~ dgamma(0.001, 0.001)
  beta_stand_dbh ~ dnorm(0, 0.001)
  beta_cdens ~ dnorm(0, 0.001)
  beta_udens ~ dnorm(0, 0.001)
  
  sigma_block ~ dgamma(0.001, 0.001)
  tau_block <- 1/sigma_block^2

  ## Predictions:
  for(m in 1:length(ud_pred)){
    
    log(out[m]) <- alpha_plot_mean + beta_udens*ud_pred[m]
    
  }
  
}

