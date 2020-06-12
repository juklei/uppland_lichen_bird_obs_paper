## bird bpo model with a variance covariance matrix for p_occ and p_det
##
## First edit: 20190326
## Last edit: 20190326
##
## Author: Julian Klein

model{
  
  ## Likelihood:
  
  ## Observational model:
  for(k in 1:nspecies){
    logit(p_det[k]) <- lp[k]
    for(i in 1:nsites){
      nseen[i,k] ~ dbin(occ_true[i,k]*p_det[k], nvisits[i,k])
      nseen_sim[i,k] ~ dbin(occ_true[i,k]*p_det[k], nvisits[i,k])
    }
  }
  
  ## Ecological process model:
  for(k in 1:nspecies){
    logit(p_occ[k]) <- lpsi[k]
    for(i in 1:nsites){
      occ_true[i,k] ~ dbern(p_occ[k])
    }
  }
  
  ## Priors:

  for(k in 1:nspecies){
    lpsi[k] <- eta[k,1]
    lp[k] <- eta[k,2]
    eta[k, 1:2] ~ dmnorm(mu_eta[], Omega[,])
  }

  for(v in 1:2){
    mu_eta[v] <- log(probs[v]/(1-probs[v]))
    probs[v] ~ dunif(0,1)
  }

  Omega[1:2, 1:2] ~ dwish(R[,], df)
  Sigma[1:2, 1:2] <- inverse(Omega[,])

  ## Predictions:
  
  ## Correlation coefficient
  rho <- Sigma[1,2]/sqrt(Sigma[1,1]* Sigma[2,2])
  
  ## Species richness:
  for(i in 1:nsites){
    richness[i] ~ sum(occ_true[i,])
  }
  scaled_rb <- (richness - mean(richness))/sd(richness)
  
}


