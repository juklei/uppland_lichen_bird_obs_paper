## model bird richness in a hierarchical model 
## 
## First edit: 20190307
## Last edit: 
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)
library(magrittr)
library(reshape2)

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

## Reduce forest variables to unique values at plot level:
pld <- unique(bpo[, c(2,9:19)])

## Create data arrays:

nvisits <- acast(bpo[, c("plot", "species", "n_visits", "obs_year")],
                 formula = plot ~ species ~ obs_year, 
                 value.var = "n_visits")
nvisits[is.na(nvisits)] <- 0

nseen <- acast(bpo[, c("plot", "species", "n_obs", "obs_year")],
               formula = plot ~ species ~ obs_year, 
               value.var = "n_obs")

## Create model data set:
data <- list(nyears = length(unique(bpo$obs_year)),
             nsites = nlevels(bpo$plot),
             nspecies = nlevels(bpo$species),
             nvisits = nvisits,
             nseen = nseen,
             canopy_density = scale(pld$PercentAbove5m),
             understory_density = scale(pld$PercentBelow5m))

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$understory_density),
                    max(data$understory_density),
                    0.05)

## Canopy density:
data$cd_pred <- seq(min(data$canopy_density), max(data$canopy_density), 0.05)

str(data)

## Inits for occ_true
T1 <- data$nseen
T1 <- ifelse(T1 > 0, 1, 0)

inits <-  list(list(occ_true = T1,
                    p_det = rep(0.4, data$nspecies),
                    alpha_mean = 2,
                    beta_ud = rep(0, data$nspecies),
                    sigma_year = 0.2,
                    sigma_spec = 2,
                    year_effect = rep(0, data$nyear),
                    spec_effect = rep(0, data$nspecies))
               )

model <- "scripts/JAGS/bird_JAGS_bpo.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  20000

update(jm, n.iter = burn.in) 

samples <- 10000
n.thin <- 5

zc <- coda.samples(jm,
                   variable.names = c("alpha_mean",
                                      "sigma_year",
                                      "sigma_spec",
                                      "beta_ud",
                                      "occ_true",
                                      "p_det"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/bpo/parameters_bird.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bird.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/bpo/diagnostics_bird.txt")

## Produce validation metrics: 
zj_val <- jags.samples(jm, 
                       variable.names = c("mean_nseen", 
                                          "mean_nseen_sim",
                                          "p_mean", 
                                          "cv_nseen", 
                                          "cv_nseen_sim", 
                                          "p_cv", 
                                          "fit", 
                                          "fit_sim",
                                          "p_fit"), 
                       n.iter = samples, 
                       thin = n.thin)

## Fit of mean:
plot(zj_val$mean_nseen, 
     zj_val$mean_nseen_sim, 
     xlab = "mean real", 
     ylab = "mean simulated", 
     cex = .05)
abline(0, 1)
p <- summary(zj_val$p_mean, mean)
text(x = 0.4, y = 0.53, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Fit of variance:
plot(zj_val$cv_nseen, 
     zj_val$cv_nseen_sim, 
     xlab = "cv real", 
     ylab = "cv simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_cv, mean)
text(x = 1.5, y = 2.1, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Overall fit:
plot(zj_val$fit, 
     zj_val$fit_sim, 
     xlab = "ssq real", 
     ylab = "ssq simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_fit, mean)
text(x = 1300, y = 1300, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)


## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("richness"),
                        n.iter = samples, 
                        thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/plot_richness_ud_bird.png", 1500, 1200, "px", res = 200)

y <- summary(zj_pred$richness, quantile, c(.025,.5,.975))$stat
x = backscale(data$ud_pred, data$understory_density)

plot(x, y[2,], 
     col="blue", 
     xlab="Understory density", 
     ylab="Richness", 
     cex = 1.4, 
     typ = "l", 
     tck = 0.03, 
     bty = "l", 
     ylim = c(20, 40)) 
polygon(c(x, rev(x)), c(y[1,], rev(y[3,])), density = 19, col = "blue", angle = 45)
lines(x,y[1,], lty="dashed", col="blue")
lines(x,y[3,], lty="dashed", col="blue")

dev.off()

## 7. Export data for fancy figures --------------------------------------------

## -------------------------------END-------------------------------------------

## To add a year effect on true occupancy / p_occ should I put that on alpha
## in the explanatory model on p_occ? try it!

## Change how you model the year random effect, e.g. mean +- epsilon. 
## Think about wether this is even 
## required given that you already have the whole model structured according 
## to the years.

## Add log(obs_time) as an offset term!