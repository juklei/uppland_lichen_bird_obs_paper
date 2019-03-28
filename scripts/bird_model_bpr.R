## model bird richness in a hierarchical model 
## 
## First edit: 20190307
## Last edit: 20190314
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
dir.create("results/bpo")
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
bpo$seen <- ifelse(bpo$n_obs > 0, 1, 0)
bpo[, c("richness", "obs_time_total") := 
      list(sum(seen), sum(obs_time)*max(n_visits)/nlevels(species)), 
    by = c("plot", "obs_year")]

## Reduce to unique rows:
bpr <- unique(bpo[, c(2:4,9:19,21:22)])

## Create model data set:
data <- list(nobs = nrow(bpr),
             nsites = nlevels(bpr$plot),
             plot = as.numeric(bpr$plot),
             obs_year = as.numeric(bpr$obs_year),
             years = unique(as.numeric(bpr$obs_year)),
             richness = bpr$richness,
             obs_time = bpr$obs_time_total,
             cd = scale(log(bpr$PercentAbove5m)),
             ud = scale(log(bpr$PercentBelow5m)),
             stand_dbh = scale(bpr$average_dbh_all_alive))

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$ud), max(data$ud), 0.05)

## Canopy density:
data$cd_pred <- seq(min(data$cd), max(data$cd), 0.05)

## Stand dbh:
data$stand_dbh_pred <- seq(min(data$stand_dbh), max(data$stand_dbh), 0.5)

str(data)

inits <-  list(list(richness_true = data$richness,
                    alpha_p_det = 0.5,
                    beta_obs_time = 0.5,
                    alpha_plot_mean = 2,
                    beta_stand_dbh = 0.2,
                    beta_cd = 0.2,
                    beta_ud = 0.2,
                    sigma_year = 2,
                    sigma_plot = 2)
               )

model <- "scripts/JAGS/bird_JAGS_bpr.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 10000
n.thin <- 5

zc <- coda.samples(jm,
                   variable.names = c("alpha_p_det",
                                      "beta_obs_time",
                                      "alpha_plot_mean",
                                      "beta_stand_dbh",
                                      "beta_cd",
                                      "beta_ud",
                                      "sigma_year",
                                      "sigma_plot"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/bpo/parameters_bpr.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bpr.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/bpo/diagnostics_bpr.txt")

## Produce validation metrics: 
zj_val <- jags.samples(jm, 
                       variable.names = c("mean_richness", 
                                          "mean_richness_sim",
                                          "p_mean", 
                                          "cv_richness", 
                                          "cv_richness_sim", 
                                          "p_cv", 
                                          "fit", 
                                          "fit_sim",
                                          "p_fit"), 
                       n.iter = samples, 
                       thin = n.thin)

## Fit of mean:
plot(zj_val$mean_richness, 
     zj_val$mean_richness_sim, 
     xlab = "mean real", 
     ylab = "mean simulated", 
     cex = .05)
abline(0, 1)
p <- summary(zj_val$p_mean, mean)
text(x = 0.5, y = 0.66, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Fit of variance:
plot(zj_val$cv_richness, 
     zj_val$cv_richness_sim, 
     xlab = "cv real", 
     ylab = "cv simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_cv, mean)
text(x = 1.5, y = 1.9, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Overall fit:
plot(zj_val$fit, 
     zj_val$fit_sim, 
     xlab = "ssq real", 
     ylab = "ssq simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_fit, mean)
text(x = 1300, y = 1200, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("r_ud", 
                                           "r_cd", 
                                           "r_stand_dbh"),
                        n.iter = samples, 
                        thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/bpr_ud.png", 1500, 1200, "px", res = 200)

y <- summary(zj_pred$r_ud, quantile, c(.025,.5,.975))$stat
x = exp(backscale(data$ud_pred, data$ud))

plot(x, y[2,], 
     col="blue", 
     xlab="Understory density", 
     ylab="Richness", 
     cex = 1.4, 
     typ = "l", 
     tck = 0.03, 
     bty = "l", 
     ylim = c(5, 20)) 
polygon(c(x, rev(x)), c(y[1,], rev(y[3,])), density = 19, col = "blue", angle = 45)
lines(x,y[1,], lty="dashed", col="blue")
lines(x,y[3,], lty="dashed", col="blue")

dev.off()

## 7. Export data for fancy figures --------------------------------------------

## -------------------------------END-------------------------------------------
