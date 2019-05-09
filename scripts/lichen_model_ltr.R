## model lichen richness in a hierarchical model with stem diameter and tree
## species at the tree level and vegetation density and stand age at the plot
## level. block is the grouping variable
## 
##
## First edit: 20190125
## Last edit: 20190426
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
             pine = ifelse(ltr$tree_sp == "Ps", 1, 0),
             spruce = ifelse(ltr$tree_sp == "Pa", 1, 0),
             stem_dbh = scale(ltr$tree_dbh),
             stand_dbh = scale(plu$average_dbh_all_alive),
             cd = scale(plu$PercentAbove5m),
             ud = scale(plu$PercentBelow5m),
             mu_p = 0.95)

## Check for explanatory variable correlation:

T1 <- plu
T1$plot <- as.numeric(T1$plot)

T2 <- data.frame("dec" = ifelse(data$pine == 0 & data$spruce == 0, 1, 0),
                 data[c(2, 5:7)]) 

cor(merge(T2, T1, by = "plot"))[-1, 5:8]

## Add prediction data:

## Stem dbh:
data$stem_dbh_pred <- seq(min(data$stem_dbh), max(data$stem_dbh), 0.05)

## Understorey density:
data$ud_pred <- seq(min(data$ud), max(data$ud), 0.05)

## Canopy density:
data$cd_pred <- seq(min(data$cd), max(data$cd), 0.05)

## Canopy density:
data$stand_dbh_pred <- seq(min(data$stand_dbh), max(data$stand_dbh), 0.05)

str(data)

inits <- list(list(p = rep(0.8, data$nobs),
                   richness_true = rep(8, data$nobs),
                   # alpha_p = -1,
                   # beta_p = 0.1,
                   sigma_p = 0.05,
                   beta_pine = 0.5,
                   beta_spruce = 0.5,
                   beta_stem_dbh = 0.5,
                   alpha_plot = rep(2, data$nplot),
                   sigma_plot = 2,
                   alpha_plot_mean = 1,
                   beta_stand_dbh = 0.2,
                   beta_cd = -0.2,
                   beta_ud = -0.2),
              list(p = rep(0.8, data$nobs),
                   richness_true = rep(8, data$nobs),
                   # alpha_p = -1,
                   # beta_p = 0.1,
                   sigma_p = 0.05,
                   beta_pine = -0.5,
                   beta_spruce = -0.5,
                   beta_stem_dbh = -0.5,
                   alpha_plot = rep(3, data$nplot),
                   sigma_plot = 1,
                   alpha_plot_mean = 10,
                   beta_stand_dbh = -0.2,
                   beta_cd = 0.2,
                   beta_ud = 0.2),
              list(p = rep(0.8, data$nobs),
                   richness_true = rep(9, data$nobs),
                   # alpha_p = -1,
                   # beta_p = 0.1,
                   sigma_p = 0.05,
                   beta_pine = 0,
                   beta_spruce = 0,
                   beta_stem_dbh = 0,
                   alpha_plot = rep(1, data$nplot),
                   sigma_plot = 1,
                   alpha_plot_mean = 5,
                   beta_stand_dbh = 0,
                   beta_cd = 0,
                   beta_ud = 0))

model <- "scripts/JAGS/lichen_JAGS_ltr.R"

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
                   variable.names = c(# "alpha_p",
                                      # "beta_p",
                                      "sigma_p",
                                      "alpha_plot_mean",
                                      "beta_pine",
                                      "beta_spruce",
                                      "beta_stem_dbh",
                                      "sigma_plot",
                                      "beta_stand_dbh",
                                      "beta_cd",
                                      "beta_ud"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lichen.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_lichen.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lichen.txt")

gelman.diag(zc)

## Useful functions:
# ecdf(zj$mean_out)(0)
# coda.matrix <- as.matrix(zc[[1]])
# head(coda.matrix)

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
# text(x = 7, y = 10.5, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Fit of variance:
# plot(zj_val$cv_richness, 
#      zj_val$cv_richness_sim, 
#      xlab = "cv real", 
#      ylab = "cv simulated", 
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_cv, mean)
# text(x = .25, y = .335, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Overall fit:
# plot(zj_val$fit, 
#      zj_val$fit_sim, 
#      xlab = "ssq real", 
#      ylab = "ssq simulated", 
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_fit, mean)
# text(x = 480, y = 650, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("r_ud",
                                           "r_cd",
                                           "r_stand_dbh",
                                           # "ud_mean", 
                                           # "cd_mean", 
                                           "spruce_mean",
                                           "pine_mean",
                                           "dec_mean"),
                        n.iter = samples, 
                        thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/rl_dbh.png", 1500, 1200, "px", res = 200)

y <- summary(zj_pred$r_ud, quantile, c(.025,.5,.975))$stat-1
x = data$ud_pred#backscale(data$ud_pred, data$ud)

plot(x, y[2,], 
     col="blue", 
     xlab="Understory density", 
     ylab="Richness per tree", 
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

export_l <- zj_pred
export_l$ud_pred <- data$ud_pred#backscale(data$ud_pred, data$ud)
export_l$cd_pred <- data$cd_pred#backscale(data$cd_pred, data$cd)
export_l$stand_dbh_pred <- data$stand_dbh_pred#backscale(data$stand_dbh_pred, data$stand_dbh)

save(export_l, file = "clean/lichen_pred.rdata")

## -------------------------------END-------------------------------------------
