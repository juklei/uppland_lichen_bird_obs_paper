## JAGS model for species_nr on nr_skarm

model{

  ## The observation:
	for(i in 1:length(species_nr)){
	
	  species_nr[i] ~ dpois(lambda[i])
		log(lambda[i]) <- alpha.mean[plot[i]] + 
		                    year.eff[year[i]] + 
		                  b1*nr_skarm_log_cent[i] +
		                  b2*dbin_35to50[i] +
		                  b3*dbin_35to50*nr_skarm_log_cent[i]
		
	}

  ## Block/plot group effect:
  for(j in 1:n_plot){
    
    alpha.mean[j] ~ dnorm(alpha.plot[block[j]], tau.alpha.plot)
  
    for(k in 1:n_block){alpha.plot[k] ~ dnorm(alpha.block, tau.alpha.block)}
  
  }
  
  ## Year group effect:
  for(l in 1:n_year){year.eff[l] ~ dnorm(0, tau.year.eff)}
  
  #priors:
  
	alpha.block ~ dnorm(1, .05)
	
	sigma.alpha.plot ~ dgamma(0.001, 0.001)
	sigma.alpha.block ~ dgamma(0.001, 0.001)
	sigma.year.eff ~ dgamma(0.001, 0.001)

	tau.alpha.plot <- 1/sigma.alpha.plot^2
	tau.alpha.block <- 1/sigma.alpha.block^2
	tau.year.eff <- 1/sigma.year.eff^2

	b1 ~ dnorm(0,.05)
	b2 ~ dnorm(0,.05)
	b3 ~ dnorm(0,.05)
	
	#prediction:
	for(m in 1:length(dbin_35to50_pred)) {
	  
	  for(n in 1:length(nr_skarm_log_cent_pred)) {
	    
	    log(output[n, m]) <- alpha.mean + 
	                         b1*nr_skarm_log_cent_pred[n] +
	                         b2*dbin_35to50_pred[m] +
	                         b3*nr_skarm_log_cent_pred[n]*dbin_35to50_pred[m]
	    
	    }
	  
	  }
	
}


## MUST have at least 1 blank line after model closes
