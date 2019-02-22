## lichen ltr model
##
## First edit: 20190129
## Last edit: 20190222
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
                      beta_pine*pine[i] + 
                      beta_spruce*spruce[i] +
                      beta_aspen*aspen[i] +
                      beta_oak*oak[i] +
                      beta_alder*alder[i] +
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
  sigma_p ~ dunif(0, 0.1)
  
  beta_pine ~ dnorm(0, 0.01)
  beta_spruce ~ dnorm(0, 0.01)
  beta_aspen ~ dnorm(0, 0.01)
  beta_oak ~ dnorm(0, 0.01)
  beta_alder ~ dnorm(0, 0.01)
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

  #R2 <- 1 - (mean(fit_sim)/fit)
  
  ## Predictions:

  ## Stem dbh:
  for(n in 1:length(stem_dbh_pred)){

    log(s_dbh_bt[n]) <- alpha_plot_mean + beta_stem_dbh*stem_dbh_pred[n]
    log(s_dbh_ps[n]) <- alpha_plot_mean + beta_stem_dbh*stem_dbh_pred[n] + 
                        beta_pine
    log(s_dbh_pa[n]) <- alpha_plot_mean + beta_stem_dbh*stem_dbh_pred[n] + 
                        beta_spruce
    log(s_dbh_pt[n]) <- alpha_plot_mean + beta_stem_dbh*stem_dbh_pred[n] + 
                        beta_aspen
    log(s_dbh_qr[n]) <- alpha_plot_mean + beta_stem_dbh*stem_dbh_pred[n] + 
                        beta_oak
    log(s_dbh_ag[n]) <- alpha_plot_mean + beta_stem_dbh*stem_dbh_pred[n] + 
                        beta_alder
    s_dbh_mean[n] <- (s_dbh_bt[n] + 
                      s_dbh_ps[n] + 
                      s_dbh_pa[n] + 
                      s_dbh_pt[n] + 
                      s_dbh_qr[n] + 
                      s_dbh_ag[n])/6

  }

  ## Understorey density:
  for(o in 1:length(ud_pred)){

    log(ud_bt[o]) <- alpha_plot_mean + beta_udens*ud_pred[o]
    log(ud_ps[o]) <- alpha_plot_mean + beta_udens*ud_pred[o] + beta_pine
    log(ud_pa[o]) <- alpha_plot_mean + beta_udens*ud_pred[o] + beta_spruce
    log(ud_pt[o]) <- alpha_plot_mean + beta_udens*ud_pred[o] + beta_aspen
    log(ud_qr[o]) <- alpha_plot_mean + beta_udens*ud_pred[o] + beta_oak
    log(ud_ag[o]) <- alpha_plot_mean + beta_udens*ud_pred[o] + beta_alder
    ud_mean[o] <- (ud_bt[o] + 
                   ud_ps[o] + 
                   ud_pa[o] + 
                   ud_pt[o] +
                   ud_qr[o] + 
                   ud_ag[o])/6

  }

  ## Canopy density:
  for(p in 1:length(cd_pred)){

    log(cd_bt[p]) <- alpha_plot_mean + beta_cdens*cd_pred[p]
    log(cd_ps[p]) <- alpha_plot_mean + beta_cdens*cd_pred[p] + beta_pine
    log(cd_pa[p]) <- alpha_plot_mean + beta_cdens*cd_pred[p] + beta_spruce
    log(cd_pt[p]) <- alpha_plot_mean + beta_cdens*cd_pred[p] + beta_aspen
    log(cd_qr[p]) <- alpha_plot_mean + beta_cdens*cd_pred[p] + beta_oak
    log(cd_ag[p]) <- alpha_plot_mean + beta_cdens*cd_pred[p] + beta_alder
    cd_mean[p] <- (cd_bt[p] + 
                   cd_ps[p] + 
                   cd_pa[p] + 
                   cd_pt[p] +
                   cd_qr[p] + 
                   cd_ag[p])/6

  }

  ## Tree species:
  log(birch_mean) <- alpha_plot_mean
  log(pine_mean) <- alpha_plot_mean + beta_pine
  log(spruce_mean) <- alpha_plot_mean + beta_spruce
  log(aspen_mean) <- alpha_plot_mean + beta_aspen
  log(oak_mean) <- alpha_plot_mean + beta_oak
  log(alder_mean) <- alpha_plot_mean + beta_alder
  
}


