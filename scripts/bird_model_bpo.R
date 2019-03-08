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
library(reshape)

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

## Reduce forestry variables to unique values at plot level:
pld <- unique(bpo[, c(2,9:19)])

## Try everything with data for one year first:
bpo <- bpo[bpo$obs_year == 2017, ]

## Create model data set:
data <- list(nsites = nlevels(bpo$plot),
             nspecies = nlevels(bpo$species),
             nvisits = as.matrix(cast(bpo[, c("plot", "species", "n_visits")], 
                                      formula = plot ~ species, 
                                      value = "n_visits")),
             nseen = as.matrix(cast(bpo[, c("plot", "species", "n_obs")], 
                                    formula = plot ~ species, 
                                    value = "n_obs")),
             canopy_density = scale(pld$PercentAbove5m),
             understory_density = scale(pld$PercentBelow5m)
             )

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$understory_density),
                    max(data$understory_density),
                    0.05)

## Canopy density:
data$cd_pred <- seq(min(data$canopy_density),  max(data$canopy_density), 0.05)

str(data)

T1 <- data$nseen
T1 <- ifelse(T1 > 0, 1, 0)
inits <-  list(list(occ_true = T1,
                    p_det = rep(0.4, data$nspecies),
                    p_occ = rep(0.4, data$nspecies)))

model <- "scripts/JAGS/bird_JAGS_bpo.R"

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
                   variable.names = c("occ_true",
                                      "p_det",
                                      "p_occ"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/bpo/parameters_bird.txt")

zj_results <- jags.samples(jm, 
                           variable.names = c("p_det"),
                           n.iter = samples, 
                           thin = n.thin)

T2 <- summary(zj_results$p_det, quantile, c(.025,.5,.975))$stat
colnames(T2) <- unique(bpo$species)
write.csv(t(T2), "results/bpo/detection_probability.csv")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bird.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/bpo/diagnostics_bird.txt")

gelman.diag(zc)

## Useful functions:
# ecdf(zj$mean_out)(0)
# coda.matrix <- as.matrix(zc[[1]])
# head(coda.matrix)

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
                                          "p_fit",
                                          "R2"), 
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
text(x = 7, y = 10.5, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Fit of variance:
plot(zj_val$cv_richness, 
     zj_val$cv_richness_sim, 
     xlab = "cv real", 
     ylab = "cv simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_cv, mean)
text(x = .25, y = .335, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Overall fit:
plot(zj_val$fit, 
     zj_val$fit_sim, 
     xlab = "ssq real", 
     ylab = "ssq simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_fit, mean)
text(x = 480, y = 650, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("richness"),
                        n.iter = samples, 
                        thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/plot_richness_cd.png", 1500, 1200, "px", res = 200)

y <- summary(zj_pred$cd_mean, quantile, c(.025,.5,.975))$stat
x = backscale(data$cd_pred, data$canopy_density)

plot(x, y[2,], 
     col="blue", 
     xlab="Canopy density", 
     ylab="Richness per averge tree", 
     cex = 1.4, 
     typ = "l", 
     tck = 0.03, 
     bty = "l", 
     ylim = c(5, 15)) 
polygon(c(x, rev(x)), c(y[1,], rev(y[3,])), density = 19, col = "blue", angle = 45)
lines(x,y[1,], lty="dashed", col="blue")
lines(x,y[3,], lty="dashed", col="blue")

dev.off()

## 7. Export data for fancy figures --------------------------------------------

export_l <- zj_pred
export_l$stem_dbh_pred <- backscale(data$stem_dbh_pred, data$stem_dbh)
export_l$ud_pred <- backscale(data$ud_pred, data$understory_density)
export_l$cd_pred <- backscale(data$cd_pred, data$canopy_density)

save(export_l, file = "clean/bird_pred.rdata")

## -------------------------------END-------------------------------------------
