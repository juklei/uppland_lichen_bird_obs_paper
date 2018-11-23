## JAGS model for species_nr on nr_skarm

model{

  ## The observation:
	for(i in 1:length(species_nr)){
	
	  species_nr[i] ~ dpois(lambda[i])
		log(lambda[i]) <- alpha.mean + 
		                    plot.eff[plot[i]] + 
		                    block.eff[block[i]] + 
		                    year.eff[obs_year[i]] + 
		                    obs.eff[observer[i]] + 
		                  b1*nr_skarm[i]
		
	}
	
  ## Plot group effect:
  for(j in 1:n_plot){plot.eff[j] ~ dnorm(0, tau.plot.eff)}
  
  ## Block group effect:
  for(k in 1:n_block){block.eff[k] ~ dnorm(0, tau.block.eff)}
  
  ## Year group effect:
  for(l in 1:n_year){year.eff[l] ~ dnorm(0, tau.year.eff)}
  
  ## Observer group effect:
  for(m in 1:n_observer){obs.eff[m] ~ dnorm(0, tau.obs.eff)}
  
  #priors:
  
	alpha.mean ~ dnorm(1, .05)
	
	sigma.plot.eff ~ dgamma(0.001, 0.001)
	sigma.block.eff ~ dgamma(0.001, 0.001)
	sigma.year.eff ~ dgamma(0.001, 0.001)
	sigma.obs.eff ~ dgamma(0.001, 0.001)
	
	tau.plot.eff <- 1/sigma.plot.eff^2
	tau.block.eff <- 1/sigma.block.eff^2
	tau.year.eff <- 1/sigma.year.eff^2
	tau.obs.eff <- 1/sigma.obs.eff^2

	b1 ~ dnorm(0,.05)
	
	#prediction:
	for(ii in 1:126) {log(output[ii]) <- alpha.mean + b1*(ii-1)}
	
}


## MUST have at least 1 blank line after model closes
