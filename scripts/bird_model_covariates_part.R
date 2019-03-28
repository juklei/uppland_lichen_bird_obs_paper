## bird and lichen richness modeled simultaniously with covariates
## 
## First edit: 20190326
## Last edit: 20190326
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
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

bpo <- read.csv("clean/bpo_50_35m.csv")

dir("clean")

load("clean/rb_2016_all_35m.rdata")
bpr_2016 <- zj_pred
load("clean/rb_2017_all_35m.rdata")
bpr_2017 <- zj_pred
load("clean/rb_2018_all_35m.rdata")
bpr_2018 <- zj_pred

## 4. The model ----------------------------------------------------------------

## Create data set with all needed vars:

d_2016 <- data.frame(r_mean = summary(bpr_2016$richness, mean)$stat,
                     r_sd = summary(bpr_2016$richness, sd)$stat,
                     obs_year = 2016,
                     plot = bpr_2016$plotnames)

d_2017 <- data.frame(r_mean = summary(bpr_2017$richness, mean)$stat,
                     r_sd = summary(bpr_2017$richness, sd)$stat,
                     obs_year = 2017,
                     plot = bpr_2017$plotnames)

d_2018 <- data.frame(r_mean = summary(bpr_2018$richness, mean)$stat,
                     r_sd = summary(bpr_2018$richness, sd)$stat,
                     obs_year = 2018,
                     plot = bpr_2018$plotnames)

d_all <- rbind(d_2016, d_2017, d_2018)

d_all <- merge(d_all, unique(bpo[, c(2,9:19)]), all.x = TRUE, by = "plot")

## Create model data set:
data <- list(nobs = nrow(d_all),
             nsites = nlevels(d_all$plot),
             sites = as.numeric(d_all$plot),
             obs_year = as.numeric(d_all$obs_year),
             years = unique(d_all$obs_year),
             r_mean = d_all$r_mean,
             r_sd = d_all$r_sd,
             cd = scale(d_all$PercentAbove5m),
             ud = scale(d_all$PercentBelow5m),
             nr_gran = scale(d_all$nr_gran),
             stand_dbh = scale(d_all$average_dbh_all_alive))

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$ud), max(data$ud), 0.05)

## Canopy density:
data$cd_pred <- seq(min(data$cd), max(data$cd), 0.05)

## Stand dbh:
data$stand_dbh_pred <- seq(min(data$stand_dbh), max(data$stand_dbh), 0.05)

str(data)

inits <-  list(list(alpha = 15,
                    beta_ud = 0.5,
                    beta_gran = 0.5,
                    sigma_year = 0.5,
                    sigma_site = 0.5)
               )

model <- "scripts/JAGS/bird&lichen_combined_JAGS_pr.R"

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
                   variable.names = c("beta_ud",
                                      "beta_gran",
                                      "alpha",
                                      "richness",
                                      "tau_site",
                                      "tau_year"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_combined.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_combined.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_combined.txt")

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("r_ud"),
                        n.iter = samples, 
                        thin = n.thin)

## Plotting prediction & 95% CIs using polygon:

png("figures/plot_richness_ud_combined.png", 1500, 1200, "px", res = 200)

y <- summary(zj_pred$r_ud, quantile, c(.025,.5,.975))$stat
x = backscale(data$ud_pred, data$ud)

plot(x, y[2,], 
     col="blue", 
     xlab="Understory density", 
     ylab="Richness", 
     cex = 1.4, 
     typ = "l", 
     tck = 0.03, 
     bty = "l", 
     ylim = c(0, 30)) 
polygon(c(x, rev(x)), c(y[1,], rev(y[3,])), density = 19, col = "blue", angle = 45)
lines(x,y[1,], lty="dashed", col="blue")
lines(x,y[3,], lty="dashed", col="blue")

dev.off()

## 7. Export data for fancy figures --------------------------------------------

## -------------------------------END-------------------------------------------
