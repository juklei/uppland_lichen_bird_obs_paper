## summed bird and lichen richness
## 
## First edit: 20190716
## Last edit: 20200612
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

forest <- read.csv("data/forest_data.csv")

dir("clean")

load("clean/rb_2017.rdata")
bpr_2017 <- zj_bird
load("clean/rb_2018.rdata")
bpr_2018 <- zj_bird
load("clean/rl.rdata")

## 4. The model ----------------------------------------------------------------

## Create data set with all needed vars:

l_2018 <- data.table(r_mean_sc = apply(zj_lichen$scaled_rl, 2, mean),
                     r_var_sc = apply(zj_lichen$scaled_rl, 2, sd)^2,
                     plot = zj_lichen$plotnames)

## Chose the following code section for the 2017 bird data:
b_2017 <- data.table(r_mean_sc = summary(bpr_2017$scaled_rb, mean)$stat,
                     r_var_sc = summary(bpr_2017$scaled_rb, sd)$stat^2,
                     plot = bpr_2017$plotnames)
d_summed <- rbind(l_2018, b_2017)
d_summed <- d_summed[, list(r_mean = sum(r_mean_sc), r_sd = sqrt(sum(r_var_sc))),
                     by = "plot"]
## --

## Chose the following code section for the 2017 bird data:
# b_2018 <- data.frame(r_mean_sc = summary(bpr_2018$scaled_rb, mean)$stat,
#                      r_var_sc = summary(bpr_2018$scaled_rb, sd)$stat^2,
#                      plot = bpr_2018$plotnames)
# l_2018 <- droplevels(l_2018[l_2018$plot %in% b_2018$plot, ])
# d_summed <- rbind(b_2018, l_2018)
# d_summed <- d_summed[, list(r_mean = sum(r_mean_sc), r_sd = sqrt(sum(r_var_sc))),
#                      by = "plot"]
##--

## Combine with forest data:
d_summed <- merge(d_summed, forest, all.x = TRUE, by = "plot")

## Create model data set:
data <- list(nobs = nrow(d_summed),
             r_mean = d_summed$r_mean,
             r_sd = d_summed$r_sd,
             od = scale(d_summed$od3), ## Chose height level break !!!!!!!!!!!!!
             ud = scale(d_summed$ud3), ## Chose height level break !!!!!!!!!!!!!
             dbh = scale(d_summed$DBH_ground))

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$ud), max(data$ud), 0.05)

## Overstory density:
data$od_pred <- seq(min(data$od), max(data$od), 0.05)

str(data)

inits <-  list(list(alpha = 5,
                    beta_ud = 0.5, beta_od = 0.5, beta_dbh = 0.5,
                    plot_sd = 0.8),
               list(alpha = -2,
                    beta_ud = 0.7, beta_od = 0.1, beta_dbh = 0,
                    plot_sd = 0.1),
               list(alpha = 0,
                    beta_ud = 1, beta_od = 1, beta_dbh = -0.5,
                    plot_sd = 1.1))

model <- "scripts/JAGS/lichen&bird_summed.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 10000, 
                 inits = inits, 
                 n.chains = 3) 

burn.in <-  90000

update(jm, n.iter = burn.in) 

samples <- 100000
n.thin <- 50

zc <- coda.samples(jm,
                   variable.names = c("alpha",
                                      "beta_ud",
                                      "beta_od",
                                      "beta_dbh",
                                      "plot_sd"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_summed_2017_3m.txt") ## Heightbreak name !!!!!!!!

## 5. Validate the model and export validation data and figures ----------------

## Look at correlations of parameters:
zc_mat <- as.matrix(zc[[1]]); dimnames(zc_mat)
plot(zc_mat[, 2], zc_mat[, 3])

pdf("figures/plot_summed_2017_3m.pdf") ## Heightbreak name !!!!!!!!!!!!!!!!!!!!!
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_summed_2017_3m.txt") ## Heightbreak name !!!!!!!

# Produce validation metrics:
zj_val <- jags.samples(jm,
                       variable.names = c("mean_r_mean", "mean_r_sim", "p_mean",
                                          "sd_r_mean", "sd_r_sim", "p_sd", 
                                          "fit", "fit_sim", "p_fit"),
                       n.iter = samples,
                       thin = n.thin)

## Fit of mean:
plot(zj_val$mean_r_mean,
     zj_val$mean_r_sim,
     xlab = "mean real",
     ylab = "mean simulated",
     cex = .05)
abline(0, 1)
mean(zj_val$p_mean)

## Fit of variance:
plot(zj_val$sd_r_mean,
     zj_val$sd_r_sim,
     xlab = "sd real",
     ylab = "sd simulated",
     cex = .05)
abline(0,1)
mean(zj_val$p_sd)

## Overall fit:
plot(zj_val$fit,
     zj_val$fit_sim,
     xlab = "ssq real",
     ylab = "ssq simulated",
     cex = .05)
abline(0,1)
mean(zj_val$p_fit)

## 6. Produce and export data for fancy figures --------------------------------

zj_pred <- jags.samples(jm,
                        variable.names = c("r_ud", "r_od"),
                        n.iter = samples,
                        thin = n.thin)

zj_pred_2017_sum <- zj_pred
zj_pred_2017_sum$ud <- backscale(data$ud_pred, data$ud)
zj_pred_2017_sum$od <- backscale(data$od_pred, data$od)
save(zj_pred_2017_sum, 
     file = "clean/summed_pred_2017_3m.rdata") ## Heightbreak name !!!!!!!!!!!!!

## -------------------------------END-------------------------------------------
