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
data$cd_pred <- seq(min(data$canopy_density),  max(data$canopy_density), 0.05)

str(data)

## Inits for occ_true
T1 <- data$nseen
T1 <- ifelse(T1 > 0, 1, 0)

## Inits for p_occ
T2 <- as.matrix(data$nseen[1,,])
T2[] <- 0.4 

inits <-  list(list(occ_true = T1,
                    p_det = rep(0.4, data$nspecies),
                    p_occ = T2,
                    alpha = rep(0, data$nspecies),
                    beta = rep(0, data$nspecies))
               )

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

T3 <- summary(zj_results$p_det, quantile, c(.025,.5,.975))$stat
colnames(T3) <- unique(bpo$species)
write.csv(t(T3), "results/bpo/detection_probability.csv")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bird.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/bpo/diagnostics_bird.txt")

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("richness"),
                        n.iter = samples, 
                        thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/plot_richness_ud.png", 1500, 1200, "px", res = 200)

y <- summary(zj_pred$ud_mean, quantile, c(.025,.5,.975))$stat
x = backscale(data$ud_pred, data$canopy_density)

plot(x, y[2,], 
     col="blue", 
     xlab="Understory density", 
     ylab="Richness", 
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

## -------------------------------END-------------------------------------------

## Notes:

# Why is p_occ not also indexed by site? Do I have to index it by year?

# Adding explanatory variables on p_occ does not work. Do I need to fully index
# also the explanatory data, e.g. understorey density[i,k,y]?

# The nesting of year, site, species; is it correct?

# Even though nseen for some plots in 2018 are NA, it estimates true occupancy
# also for those plots. I am assuming its the priors? Is the solution to set the
# priors to 0? Why is it not enough to state that the nvisits for that spot is
# 0 to not get an estimate for true occupancy?

## To add a year effect on true occupancy / p_occ should I put that on alpha
## in the explanatory model on p_occ?

