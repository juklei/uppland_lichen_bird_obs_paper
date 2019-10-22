## model lichen richness per stand with species accumulation curves
## The species accumulation curve is a michaelis-menten saturation curve
## We extract the asymptote richness and model it together with bird richness
##
## First edit: 20190605
## Last edit: 20190617
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
library(gmodels)

## 2. Define or source functions used in this script ---------------------------

dir.create("results")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## 3. Load and explore data ----------------------------------------------------

dir("clean")

load("clean/species_accumulation_data.rda")
load("temp/plot_order.rda")
str(sad)

## 4. The model ----------------------------------------------------------------

## Calculate the number of trees by plot:
ntree <- apply(sad, 3, function(x) sum(!is.na(x[,2])))

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nplot = dim(sad)[3],
             ntree = ntree,
             obs = sad)

str(data)

## Prepare inits:
plot_richness <- sad[1,,]
plot_richness[] <- 20
sat_speed <- sad[1,,]
sat_speed[] <- 5

inits <- list(list(plot_richness = plot_richness,
                   lambda_rich = rep(10, data$nrep),
                   lambda_sat = rep(5, data$nrep),
                   sat_speed = sat_speed
                   ))

model <- "scripts/JAGS/lichen_JAGS_lpsac_part.R"

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
                   variable.names = "plot_richness",
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lpsac_part.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_lpsac_part.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lpsac_part.txt")

## Produce validation metrics:
zj_val <- jags.samples(jm,
                       variable.names = "obs_pred",
                       n.iter = 1000,
                       thin = 10)

pred <- summary(zj_val$obs_pred, quantile, c(.025,.5,.975))$stat
x = 0:50

median_mean <- apply(pred[2,,,], c(1,3), mean)
ci_max <- apply(pred[3,,,], c(1,3), max)
ci_min <- apply(pred[1,,,], c(1,3), min)

dev.off()

pdf("figures/sim_vs_obs.pdf")

par(mfrow = c(3, 2))

for(i in 1:data$nplot) {

  plot(x, c(0, ci_max[,i]), 
       lty = "dashed", 
       col = "blue", 
       xlab = "tree nr", 
       ylab = "richness", 
       typ = "l")
  lines(x, c(0,median_mean[,i]), col = "blue")
  lines(x, c(0, ci_min[,i]), lty = "dashed", col = "blue")

  ## Real data:
  points(rep(which(!is.na(sad[,1,i])), dim(sad)[2]),
         na.omit(as.vector(sad[,,i])))
  
}

dev.off()

## 6. Export data --------------------------------------------------------------

zj_lichen <- jags.samples(jm,
                          variable.names = c("plot_richness", "scaled_rl"),
                          n.iter = samples,
                          thin = n.thin)

zj_lichen$plotnames <- plot_order

save(zj_lichen, file = "clean/lichen_richness.rdata")

## -------------------------------END-------------------------------------------
