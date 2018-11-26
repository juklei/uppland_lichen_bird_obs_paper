## IIn this script bird obs data used in the bird & lichen paper is modelled
## with different forest variables. A bayesian approach is used. Different 
## responses are modelled from species to community level.
##
## First edit: 20181122
## Last edit: 20181123
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(dplyr)
library(boot)
library(rjags)
library(coda)
library(corrgram)

## 2. Define or source functions used in this script ---------------------------

#source()

## 3. Load and explore data ----------------------------------------------------

dir("data")
f_plot <- read.csv("data/forest_data_uppland_plot.csv")
b_obs <- read.csv("data/occ_2016to2018.csv")
excl <- read.csv("data/non_experiment_include.csv")

head(f_plot); head(b_obs); head(excl)

## 4. Merge bird with forest data and exclude post treatment observations ------ 

bf <- merge(b_obs, f_plot, all.x = TRUE, by = c("plot", "block"))

bf <- merge(bf, excl, by = "plot")
bf <- bf[!(bf$obs_year == 2018 & bf$include_2018 == "N"), -length(bf)]

## 5. Explore and transform the data set for statistical analysis --------------

## Look at explanatory variable correlations and export:

dir.create("results")

cor(unique(na.omit(bf[, c(17:19, 22:31)]))) %>% 
  capture.output() %>% write(., "results/forest_correlations.txt")

dir.create("figures")

png("figures/forest_correlations.png", 3000, 1000, "px")

corrgram(unique(na.omit(bf[, c(17:19, 22:31)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

## Exclude predators and siskins from obs
bf <- bf[!bf$species %in% c("grona", "bergk", "ekore", "mard", "gravg"), ]

## Add distance categories:

bf$dist_bin[bf$dist <= (50/sqrt(2))] <- "0to35m"
bf$dist_bin[bf$dist > (50/sqrt(2))] <- "35to50m"
# 
# T1 <- bf
# T1$dist_bin <- "0to50"
# 
# bf <- rbind(bf, T1)

## For now, take only obstime <=15min to compare years:
## Later all observations might enter analysis by means if rarefaction curves
bf <- bf[bf$minutes_to_obs < 16, ]

## Add species numbers;

bf <- as.data.table(bf)

## per year/plot/visit;
bf[, "species_nr" := length(unique(species)), 
   by = c("plot", "visit", "obs_year", "dist_bin")] 

## and per year/plot/visit/functional group;
# bf[, "species_nr" := length(unique(species)), 
#    by = c("plot", "visit", "obs_year", "dist_bin", "func_group")] 

## 5. Reduce data set and make model species_nr or abundances ------------------

## 5.a Predict nr_species with nr. umbrella spruce -----------------------------

# Reduce data set and store as list:
red1 <- na.omit(unique(bf[, c("plot",
                              "block", 
                              "visit",
                              "observer",
                              "obs_year", 
                              "dist_bin", 
                              "species_nr", 
                              "nr_skarm"), ]))

str(red1)
hist(red1$species_nr)

## Transform data set to list for JAGS:

D1 <- as.list(NULL)

## Observation:
D1$id <- 1:nrow(red1)

## Response:
D1$species_nr <- red1$species_nr

## Fixed:
D1$nr_skarm_log_cent <- log(red1$nr_skarm+1)-mean(log(red1$nr_skarm+1))
D1$dbin_35to50 <- ifelse(red1$dist_bin == "35to50", 1, 0)

## Random:
D1$plot <- as.numeric(red1$plot)
D1$block <- as.numeric(red1$block)
D1$year <- red1$obs_year-2015

##Prediction:
D1$dbin_35to50_pred <- c(0, 1)
D1$nr_skarm_log_cent_pred <- seq(0, max(log(red1$nr_skarm+1)), 0.1) - 
  mean(log(red1$nr_skarm+1))

## Add nlevels of obs and group effects:
D1$n_obs <- nrow(red1)
D1$n_plot <- max(D1$plot)
D1$n_block <- max(D1$block)
D1$n_year <- max(D1$year)

str(D1)

## Define initial values of parameters:
inits <- list(list(alpha.mean = 1, 
                   alpha.plot = 1,
                   alpha.block = 1,
                   b1 = 0.1,
                   b2 = 0,
                   b3 = 0,
                   sigma.alpha.plot = 1,
                   sigma.alpha.block = 1,
                   sigma.year.eff = 1))

## Load the JAGS model file:
bsp_nr_nr_skarm <- "scripts/JAGS/bird_forest_JAGS_bsp_nr.R"

## Compile the model:
jm <- jags.model(bsp_nr_nr_skarm, 
                 data = D1,  
                 n.adapt = 5000,
                 inits = inits, 
                 n.chains = 1) 

## Burn-inthe model with n.iter = nr. of iterations:
update(jm, n.iter = 10000) 

## Generate samples from the posterior distribution with CODA:

samples <- 10000
n.thin <- 5 ## Takes only every fifth sample

zc <- coda.samples(jm,
                   variable.names = c("alpha.mean",
                                      "alpha.plot",
                                      "alpha.block",
                                      "b1",
                                      "b2", 
                                      "b3",
                                      "sigma.alpha.plot",
                                      "sigma.alpha.block",
                                      "sigma.year.eff"), 
                   n.iter = samples, 
                   thin = n.thin) 

## Explore:
summary(zc) #will show the mean estimate and SE and 95% CIs
plot(zc) #look at the chains to see stability and mixing


## Generate samples from the posterior distribution with JAGS

zj <- jags.samples(jm, 
                   variable.names = c("alpha.mean", "b1", "b2", "b3", "output"), 
                   n.iter = samples, 
                   thin = n.thin)

## Diagnostics for convergence:
heidel.diag(zc)
gelman.diag(zc) #needs at least 2 chains

## Export model results and predictions:

## Generate medians and quantiles that can be used for storing info, plotting etc.
pred <- summary(zj$output, quantile, c(.025,.5,.975))$stat

cbind("pred_nr_spec" = t(pred), "nr_skarm" = 0:125) %>% 
  write.csv(., "data/p.nr_spec_nr_skarm.csv", row.names = FALSE)

## 5.b Predict nr_species with xxx ---------------------------------------------

## ...

## -------------------------------END-------------------------------------------

## Quick graphing:

#plotting prediction & 95%CIs using polygon
plot(0:125, 
     pred[2,], 
     col = "blue", 
     xlab = "nr umbrella spruce", 
     ylab = "nr species", 
     cex = 1.4, 
     typ = "l", 
     tck = 0.03, 
     bty = "l")
polygon(c(0:125, rev(0:125)), 
        c(pred[1,], rev(pred[3,])), 
        density = 19, 
        col = "blue", 
        angle = 45)
lines(0:125,  pred[1,], lty = "dashed", col = "blue")
lines(0:125,  pred[3,], lty = "dashed", col = "blue")

