## model bird richness in a hierarchical model 
## 
## First edit: 20190307
## Last edit: 20190426
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)
library(magrittr)
library(reshape2)
library(data.table)

## 2. Define or source functions used in this script ---------------------------

dir.create("results")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

dir("clean")

bpo <- read.csv("clean/bpo_50.csv")
head(bpo)
str(bpo)

## 4. The model ----------------------------------------------------------------

## Calculate richness and total observation time per plot and year:
bpo <- as.data.table(bpo)
bpo$seen <- ifelse(bpo$n_obs > 0, 1, 0) ## At least two obs per plot and season
bpo[, c("richness", "obs_time_total") := 
      list(sum(seen), sum(obs_time)*max(n_visits)/nlevels(species)), 
    by = c("plot", "obs_year")]

## Reduce to unique rows:
bpr <- unique(bpo[, c(2:4,9:20,22:23)])

## Create model data set:
data <- list(nobs = nrow(bpr),
             nsites = nlevels(bpr$plot),
             plot = as.numeric(bpr$plot),
             obs_year = as.numeric(bpr$obs_year),
             years = unique(as.numeric(bpr$obs_year)),
             richness = bpr$richness,
             obs_time = bpo$obs_time*60,
             cd = scale(log(bpr$PercentAbove5m)),
             ud = scale(log(bpr$PercentBelow5m)),
             spruce = scale(log(bpr$nr_gran + 1)),
             pine = scale(log(bpr$nr_tall + 1)),
             dec = scale(log(bpr$nr_lov + 1)),
             umbrella = scale(log(bpr$nr_skarm + 1)),
             stand_dbh = scale(bpr$average_dbh_all_alive))

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$ud), max(data$ud), 0.05)

## Canopy density:
data$cd_pred <- seq(min(data$cd), max(data$cd), 0.05)

## Stand dbh:
data$stand_dbh_pred <- seq(min(data$stand_dbh), max(data$stand_dbh), 0.05)

## Nr. deciduous:
data$dec_pred <- seq(min(data$dec), max(data$dec), 0.05)

## Nr. spruce:
data$spruce_pred <- seq(min(data$spruce), max(data$spruce), 0.05)

## Nr. pine:
data$pine_pred <- seq(min(data$pine), max(data$pine), 0.05)

## Nr. umbrella:
data$umbr_pred <- seq(min(data$umbrella), max(data$umbrella), 0.05)

str(data)

inits <-  list(list(richness_true = data$richness,
                    param_obs_time = 500,
                    alpha_plot_mean = 2,
                    beta_stand_dbh = 0.2,
                    beta_cd = 0.2,
                    beta_ud = 0.2,
                    beta_spruce = 0.2,
                    beta_pine = 0.2,
                    beta_dec = 0.2,
                    beta_umbr = 0.2,
                    sigma_year = 2,
                    sigma_sites = 2),
               list(richness_true = data$richness,
                    param_obs_time = 1000,
                    alpha_plot_mean = 5,
                    beta_stand_dbh = 0.5,
                    beta_cd = 0.5,
                    beta_ud = 0.5,
                    beta_spruce = 0.5,
                    beta_pine = 0.5,
                    beta_dec = 0.5,
                    beta_umbr = 0.5,
                    sigma_year = 5,
                    sigma_sites = 5),
               list(richness_true = data$richness,
                    param_obs_time = 100,
                    alpha_plot_mean = 0,
                    beta_stand_dbh = -0.2,
                    beta_cd = -0.2,
                    beta_ud = -0.2,
                    beta_spruce = -0.2,
                    beta_pine = -0.2,
                    beta_dec = -0.2,
                    beta_umbr = -0.2,
                    sigma_year = 1,
                    sigma_sites = 1)
               )

model <- "scripts/JAGS/bird_JAGS_bpr.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 3) 

burn.in <-  100000

update(jm, n.iter = burn.in) 

samples <- 50000
n.thin <- 25

zc <- coda.samples(jm,
                   variable.names = c("param_obs_time",
                                      "alpha_plot_mean",
                                      "beta_stand_dbh",
                                      "beta_cd",
                                      "beta_ud",
                                      "beta_spruce",
                                      "beta_pine",
                                      "beta_dec",
                                      "beta_umbr",
                                      "sigma_year",
                                      "sigma_site"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_bpr_veg_age.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bpr_veg_age.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_bpr_veg_age.txt")

# ## Produce validation metrics:
# zj_val <- jags.samples(jm,
#                        variable.names = c("mean_richness",
#                                           "mean_richness_sim",
#                                           "p_mean",
#                                           "cv_richness",
#                                           "cv_richness_sim",
#                                           "p_cv",
#                                           "fit",
#                                           "fit_sim",
#                                           "p_fit"),
#                        n.iter = samples,
#                        thin = n.thin)
# 
# ## Fit of mean:
# plot(zj_val$mean_richness,
#      zj_val$mean_richness_sim,
#      xlab = "mean real",
#      ylab = "mean simulated",
#      cex = .05)
# abline(0, 1)
# p <- summary(zj_val$p_mean, mean)
# text(x = 0.5, y = 0.66, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Fit of variance:
# plot(zj_val$cv_richness,
#      zj_val$cv_richness_sim,
#      xlab = "cv real",
#      ylab = "cv simulated",
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_cv, mean)
# text(x = 1.5, y = 1.9, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Overall fit:
# plot(zj_val$fit,
#      zj_val$fit_sim,
#      xlab = "ssq real",
#      ylab = "ssq simulated",
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_fit, mean)
# text(x = 1300, y = 1200, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("r_ud",
                                           "r_cd",
                                           "r_stand_dbh",
                                           "r_dec",
                                           "r_spruce",
                                           "r_pine",
                                           "r_umbr"),
                        n.iter = samples, 
                        thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/bpr_ud_red.png", 1200, 1200, "px", res = 200)

y <- summary(zj_pred$r_ud, quantile, c(.025,.5,.975))$stat-1
#y <- y - mean(y[2,])
x = data$ud_pred#exp(backscale(data$ud_pred, data$ud))

plot(x, y[2,], 
     col="blue", 
     xlab="Understory density", 
     ylab="Richness", 
     cex = 1.4, 
     typ = "l", 
     tck = 0.03, 
     bty = "l", 
     ylim = c(-0.5, 1)) 
polygon(c(x, rev(x)), c(y[1,], rev(y[3,])), density = 19, col = "blue", angle = 45)
lines(x,y[1,], lty="dashed", col="blue")
lines(x,y[3,], lty="dashed", col="blue")

dev.off()

## 7. Export data for fancy figures --------------------------------------------

export_b_veg_age <- zj_pred
export_b_veg_age$ud_pred <- data$ud_pred#backscale(data$ud_pred, data$ud)
export_b_veg_age$cd_pred <- data$cd_pred#backscale(data$cd_pred, data$cd)
export_b_veg_age$stand_dbh_pred <- data$stand_dbh_pred#backscale(data$stand_dbh_pred, data$stand_dbh)

save(export_b_veg_age, file = "clean/bird_pred_veg_age.rdata")

# export_b_trees <- zj_pred
# export_b_trees$dec_pred <- exp(backscale(data$dec_pred, data$dec)) - 1
# export_b_trees$spruce_pred <- exp(backscale(data$spruce_pred, data$spruce)) - 1
# export_b_trees$pine_pred <- exp(backscale(data$pine_pred, data$pine) - 1)
# # export_b_trees$umbr_pred <- exp(backscale(data$umbr_pred, data$umbrella) - 1)
# 
# save(export_b_trees, file = "clean/bird_pred_trees.rdata")

## -------------------------------END-------------------------------------------
