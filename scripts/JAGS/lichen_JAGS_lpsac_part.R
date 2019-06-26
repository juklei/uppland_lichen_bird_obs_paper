## lichen lpsac model
##
## First edit: 20190605
## Last edit: 20190625
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observation model:
  for(p in 1:nplot){
    for(k in 1:nrep){
      for(i in 1:ntree[p]){
        obs[i,k,p] ~ dpois(lambda_obs[i,k,p])
        lambda_obs[i,k,p] <- (plot_richness[k,p]*i)/(sat_speed[k,p] + i)
  }}}
  
  ## Process model:
  for(p in 1:nplot){
    for(k in 1:nrep){
      plot_richness[k,p] ~ dpois(lambda_rich[k])
      sat_speed[k,p] ~ dpois(lambda_sat[k])
  }}
  
  ## Priors:
  for(k in 1:nrep){
    lambda_rich[k] ~ dgamma(0.001, 0.001)
    lambda_sat[k] ~ dgamma(0.001, 0.001)
  }
    
  ## Predictions:
  
  ## Simulation:
  for(p in 1:nplot){
    for(k in 1:nrep){
      for(m in 1:50){
        obs_pred[m,k,p] <- (plot_richness[k,p]*m)/(sat_speed[k,p] + m)
  }}}

  
  for(p in 1:nplot){
    for(k in 1:nrep){
      obs_limit[k,p] <- (plot_richness[k,p]*10E6)/(sat_speed[k,p] + 10E6)
  }}
  
  scaled_rl <- (obs_limit - mean(obs_limit))/sd(obs_limit)
  
}

