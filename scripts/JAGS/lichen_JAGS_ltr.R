## lichen ltr model
##
## First edit: 20190129
## Last edit: 20190129
##
## Author: Julian Klein

model{

## likelihood
  
## Tree level:
for(i in id){
  
  richness[i] ~ dpois(lambda[i])
  log(lambda[i]) <- alpha_plot[plot[i]] + 
                    beta_pine[i]*tree_sp_pine + 
                    beta_spruce[i]*tree_sp_spruce + 
                    beta_stem_dbh[i]*stem_dbh
  
}

## Plot level:
for(j in unique(plot)){
  
  alpha_plot[j] ~ dgamma(a_plot[j], b_plot[j])
  mu[j] <- alpha_plot_mean + block_effect[block[j]] +
           beta_stand_dbh[j]*stand_dbh +
           beta_cdens[j]*canopy_density +
           beta_udens[j]*understory_density
  
} 

## Block level:
for (k in 1:unique(block)){
  
  site.eff[k] ~ dnorm(0, tau_block)
  
}
  
#priors



#predictions

}

