## model lichen richness in a hierarchical model with stem diameter and tree
## species at the tree level and vegetation density and stand age at the plot
## level. block is the grouping variable
## 
##
## First edit: 20190125
## Last edit: 20190131
##
## Author: Julian Klein, Matt Low

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)
library(magrittr)

## 2. Define or source functions used in this script ---------------------------

dir.create("results")
dir.create("results/ltr")
dir.create("figures")

## 3. Load and explore data ----------------------------------------------------

dir("clean")

ltr <- na.omit(read.csv("clean/ltr_T_10.csv"))
head(ltr)
str(ltr)

## 4. The model ----------------------------------------------------------------

## Plot level explanatory variables need to be reduced to unique rows:
plu <- unique(ltr[, c("plot", 
                      "average_dbh_all_alive", 
                      "PercentAbove5m", 
                      "PercentBelow5m")])
nrow(plu)

## Create model data set:
data <- list(nobs = nrow(ltr),
             # block = as.numeric(ltr$block),
             # nblock = length(unique(ltr$block)),
             plot = as.numeric(ltr$plot),
             nplot = length(unique(ltr$plot)),
             richness = ltr$richness,
             tree_sp_pine = ifelse(ltr$tree_sp == "Ps", 1, 0),
             tree_sp_spruce = ifelse(ltr$tree_sp == "Pa", 1, 0),
             stem_dbh = scale(ltr$tree_dbh),
             stand_dbh = scale(plu$average_dbh_all_alive),
             canopy_density = scale(plu$PercentAbove5m),
             understory_density = scale(plu$PercentBelow5m))

## Add prediction data:
data$ud_pred <- seq(min(data$understory_density),
                    max(data$understory_density),
                    0.05)

str(data)

inits <-  list(list(p = 0.8,
                    richness_true = rep(12, data$nobs),
                    beta_pine = 0.5,
                    beta_spruce = 0.5,
                    beta_stem_dbh = 0.5,
                    alpha_plot = rep(2, data$nplot),
                    sigma_plot = 2,
                    alpha_plot_mean = 1,
                    beta_stand_dbh = 0.2,
                    beta_cdens = -0.2,
                    beta_udens = -0.2))

model <- "scripts/JAGS/lichen_JAGS_ltr.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 1000
n.thin <- 1

zc <- coda.samples(jm,
                   variable.names = c("p",
                                      "beta_pine",
                                      "beta_spruce",
                                      "beta_stem_dbh",
                                      "alpha_plot",
                                      "sigma_plot",
                                      "alpha_plot_mean",
                                      "beta_stand_dbh",
                                      "beta_cdens",
                                      "beta_udens"), 
                   n.iter = samples, 
                   thin = n.thin)


capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/ltr/parameters.txt")

plot(zc) 

## Diagnostics:

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/ltr/diagnostics.txt")

gelman.diag(zc)

## Useful functions:
# ecdf(zj$mean_out)(0)
# coda.matrix <- as.matrix(zc[[1]])
# head(coda.matrix)

## 5. Produce and export figures -----------------------------------------------

dev.off()

## Produce predictions:
zj <- jags.samples(jm, 
                   variable.names = c("out_d", "out_p", "out_s", "mean_out"), 
                   n.iter = samples, 
                   thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/plot_richness_ud.png", 1500, 1200, "px", res = 200)

pred <- summary(zj$mean_out, quantile, c(.025,.5,.975))$stat
x = data$ud_pred + mean(plu$PercentBelow5m)
y = pred
plot(x, y[2,], 
     col="blue", 
     xlab="Understorey density", 
     ylab="Richness per plot", 
     cex = 1.4, 
     typ = "l", 
     tck = 0.03, 
     bty = "l", 
     ylim = c(5, 15)) 
polygon(c(x, rev(x)), c(y[1,], rev(y[3,])), density = 19, col = "blue", angle = 45)
lines(x,y[1,], lty="dashed", col="blue")
lines(x,y[3,], lty="dashed", col="blue")

dev.off()

## -------------------------------END-------------------------------------------
