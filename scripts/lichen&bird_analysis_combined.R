## bird and lichen richness modeled simultaniously with covariates
## 
## First edit: 20190326
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

## 2. Define or source functions used in this script ---------------------------

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

forest <- read.csv("data/forest_data.csv")

## calculate per ha tree numbers and export correlations and data summary:
forest$no_dec <- forest$no_dec/(pi*.03)
forest$no_spruce <- forest$no_spruce/(pi*.03)
forest$no_pine <- forest$no_pine/(pi*.03)
forest$no_ud_spruce <- forest$no_ud_spruce/(pi*.03)
write.csv(cor(forest[, -1]), "results/correlations.csv")
capture.output(summary(forest[, -1])) %>% write(., "results/summary.txt")

## Load data from taxa specific models:
load("clean/rb_2017.rdata")
bpr_2017 <- zj_bird
load("clean/rb_2018.rdata")
bpr_2018 <- zj_bird
load("clean/rl.rdata")

## 4. The model ----------------------------------------------------------------

## Create data set with all needed vars:

l_2018 <- data.frame(r_mean_sc = apply(zj_lichen$scaled_rl, 2, mean),
                     r_sd_sc = apply(zj_lichen$scaled_rl, 2, sd),
                     obs_year = 2018,
                     organsim = "lichen",
                     plot = zj_lichen$plotnames)

## Chose the following code section for the 2017 bird data:
b_2017 <- data.frame(r_mean_sc = summary(bpr_2017$scaled_rb, mean)$stat,
                     r_sd_sc = summary(bpr_2017$scaled_rb, sd)$stat,
                     obs_year = 2017,
                     organsim = "bird",
                     plot = bpr_2017$plotnames)
d_all <- rbind(b_2017, l_2018)
## --

## Chose the following code section for the 2017 bird data:
# b_2018 <- data.frame(r_mean_sc = summary(bpr_2018$scaled_rb, mean)$stat,
#                      r_sd_sc = summary(bpr_2018$scaled_rb, sd)$stat,
#                      obs_year = 2018,
#                      organsim = "bird",
#                      plot = bpr_2018$plotnames)
# l_2018 <- droplevels(l_2018[l_2018$plot %in% b_2018$plot, ])
# d_all <- rbind(b_2018, l_2018)
## --

## Combine with forest data:
d_all <- merge(d_all, forest, all.x = TRUE, by = "plot")

## Create model data set: 
data <- list(nobs = nrow(d_all),
             org = ifelse(d_all$organsim == "lichen", 1, 0),
             r_mean = d_all$r_mean_sc,
             r_sd = d_all$r_sd_sc,
             od = scale(d_all$od3), ## Chose height level break !!!!!!!!!!!!!!!!
             ud = scale(d_all$ud3), ## Chose height level break !!!!!!!!!!!!!!!!
             dbh = scale(d_all$DBH_ground))

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$ud), max(data$ud), 0.05)

## Overstory density:
data$od_pred <- seq(min(data$od), max(data$od), 0.05)

str(data)

inits <-  list(list(alpha = 5, beta_org = 0, 
                    beta_ud = 0.5, beta_od = 0.5, beta_dbh = 0.5,
                    int_ud = 0.1, int_od = 0.2, int_dbh = 0.1,
                    plot_sd = 0.1),
               list(alpha = -2, beta_org = 0.2,
                    beta_ud = 0.7, beta_od = 0.1, beta_dbh = 0,
                    int_ud = 0, int_od = 0, int_dbh = 0,
                    plot_sd = 0.3),
               list(alpha = 0, beta_org = 0,
                    beta_ud = 1, beta_od = 1, beta_dbh = -0.5,
                    int_ud = -0.1, int_od = -0.2, int_dbh = -0.1,
                    plot_sd = 0.5))

model <- "scripts/JAGS/lichen&bird_combined.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 10000, 
                 inits = inits, 
                 n.chains = 3) 

burn.in <-  40000

update(jm, n.iter = burn.in) 

samples <- 100000
n.thin <- 50

zc <- coda.samples(jm,
                   variable.names = c("alpha", "beta_org",
                                      "beta_ud", "beta_od", "beta_dbh",
                                      "int_ud", "int_od", "int_dbh",
                                      "plot_sd",
                                      "lichen_dbh"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_combined_2017_3m_sc.txt") ## Heightbreak name !!!

## 5. Validate the model and export validation data and figures ----------------

## Look at correlations of parameters:
zc_mat <- as.matrix(zc[[1]]); dimnames(zc_mat)
plot(zc_mat[, 2], zc_mat[, 3])

pdf("figures/plot_combined_2017_3m_sc.pdf") ## Heightbreak name !!!!!!!!!!!!!!!!
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_combined_2017_3m_sc.txt") ## Heightbreak name !!

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
                        variable.names = c("rb_ud", "rl_ud", "rb_od", "rl_od"),
                        n.iter = samples,
                        thin = n.thin)
zj_pred_2017 <- zj_pred
zj_pred_2017$ud <- backscale(data$ud_pred, data$ud)
zj_pred_2017$od <- backscale(data$od_pred, data$od)
save(zj_pred_2017, 
     file = "clean/combined_pred_2017_3m_sc.rdata") ## Heightbreak name !!!!!!!!

## -------------------------------END-------------------------------------------
