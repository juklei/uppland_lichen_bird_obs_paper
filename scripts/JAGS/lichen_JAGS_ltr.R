## lichen ltr model
##
## First edit: 20190129
## Last edit: 20190129
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Tree level:
  for(i in id){
  
    richness[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha_plot[plot[i]] + 
                      beta_pine[i]*tree_sp_pine + 
                      beta_spruce[i]*tree_sp_spruce + 
                      beta_stem_dbh[i]*stem_dbh
    
  }
  
  ## Plot level:
  for(j in 1:nplot){
  
    alpha_plot[j] ~ dgamma(mu[j]^2/sigma_plot^2, mu[j]/sigma_plot^2)
    log(mu[j]) <- alpha_plot_mean + block_effect[block[j]] +
                  beta_stand_dbh[j]*stand_dbh +
                  beta_cdens[j]*canopy_density +
                  beta_udens[j]*understory_density
  
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

}

