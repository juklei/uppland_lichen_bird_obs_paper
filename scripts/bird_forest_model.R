## In this script bird obs data used in the bird & lichen paper is modelled
## with different forest variables. A bayesian approach is used. Different 
## responses are modelled from species to community level.
##
## First edit: 20181122
## Last edit: 20181123
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(reshape)
library(data.table)
library(dplyr)
library(boot)
library(rjags)
library(coda)
library(corrgram)
library(drc)
library(unmarked)

## 2. Define or source functions used in this script ---------------------------

#source()

## 3. Load and explore data ----------------------------------------------------

dir("data")
f_plot <- read.csv("data/forest_data_uppland_plot.csv")
b_obs <- read.csv("data/occ_2016to2018.csv")
excl <- read.csv("data/non_experiment_include.csv")

head(f_plot); head(b_obs); head(excl)

## Look at forest explanatory variable correlations and export:

dir.create("figures")

png("figures/forest_correlations.png", 3000, 1000, "px")

corrgram(unique(na.omit(f_plot[, c(6:8, 11:20)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

## 4. Rearrange data set so it can be used in unmarked and JAGS. --------------- 
##    Exclude post treatment data and unised species.
##    This part might later be moved to general if usage the same in several
##    analyses in this project.

## Summarise time after sunrise per visit to become an observational
## covariate in the analysis:
b_obs <- as.data.table(b_obs)
b_obs[, "mp_sunrise_visit" := mean(min_post_sunrise), 
      by = c("plot", "obs_year", "visit")]

## Make data set so we have all possible combinations of all visits
## and all species seen during the whole survey.
b_occ <- expand.grid.df(unique(b_obs[, c("block", 
                                         "plot", 
                                         "observer", 
                                         "visit", 
                                         "obs_year",
                                         "obs_time",
                                         "mp_sunrise_visit", 
                                         "dp_march")]), 
                        as.data.frame(levels(b_obs$species)))
colnames(b_occ)[9] <- "species"

## Merge with b_obs data set to aquire all information of observations
## If a speces was not observed during a visit make obs 0, otherwise NA:

## Add obs = 1 to b_obs:
b_obs$observed <- 1

## Merge all b_occ with matching b_obs:
b_occ <- merge(b_occ, b_obs, all.x = TRUE, by = c("block", 
                                                  "plot", 
                                                  "observer", 
                                                  "visit", 
                                                  "obs_year",
                                                  "obs_time",
                                                  "species",
                                                  "mp_sunrise_visit", 
                                                  "dp_march"))

## NA in observed are actually non-observations, so they will become 0:
b_occ$observed[is.na(b_occ$observed)] <- 0

## If everything went well the number of positive observations should be the 
## nrow of b_obs:
if(sum(b_occ$observed == 1) != nrow(b_obs)) print("You made a mistake!")

## Now the observations need to be put into a format that fits the unmarked 
## package and JAGS: Observational covariates are chosen here:
b_occ <- dcast(setDT(b_occ), 
               block + plot + observer + species + obs_year ~ visit, 
               value.var = c("observed", 
                             "mp_sunrise_visit",
                             "dp_march",
                             "obs_time"))

## Merge with forest data:
bf_occ <- merge(b_occ, f_plot, all.x = TRUE, by = c("plot", "block"))

## Exclude plots on which surveys were performed after treatments:
bf_occ <- merge(bf_occ, excl, by = "plot")
bf_occ <- bf_occ[!(bf_occ$obs_year == 2018 & bf_occ$include_2018 == "N"), ]

## Exclude predators and siskins from obs
bf_occ <- bf_occ[!bf_occ$species %in% c("grona", 
                                        "bergk", 
                                        "ekore", 
                                        "mard", 
                                        "gravg"), ]

## 5. Run occupancy models using unmarked and forest covariate -----------------

## Lets try it first for only one species, one observer in 2017 and nr_skarm:
rodhe <- as.data.frame(bf_occ[bf_occ$species == "rodhe" & 
                                bf_occ$obs_year == 2017, ])

## Make an unmarkedFrame:
um.rodhe <- unmarkedFrameOccu(y = as.matrix(rodhe[, c("observed_first",
                                                      "observed_second",
                                                      "observed_third",
                                                      "observed_fourth",
                                                      "observed_fifth")]),
                              siteCovs = rodhe[, c("block", "nr_skarm")],
                              obsCovs = list("date" = 
                                               rodhe[, c("dp_march_first",
                                                         "dp_march_second",
                                                         "dp_march_third",
                                                         "dp_march_fourth",
                                                         "dp_march_fifth")]))

## Fit occupancy model:

summary(um.rodhe)
occu(~ date ~ nr_skarm, data = um.rodhe)

m.um.rodhe <- occu(~ 1 ~ 1, data = um.rodhe)
backTransform(m.um.rodhe, "state")
backTransform(m.um.rodhe, "det")

## 6. Run occupancy models using JAGS and forest covariates --------------------

## Lets try it first for only one species, one observer in 2017 and nr_skarm:
rodhe <- as.data.frame(bf_occ[bf_occ$species == "rodhe" & 
                                bf_occ$obs_year == 2017, ])

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

## N of groups:
D1$n_plots <- max(D1$plot)
D1$n_years <- max(D1$year)

##Prediction:
D1$dbin_35to50_pred <- c(0, 1)
D1$nr_skarm_log_cent_pred <- seq(0, max(log(red1$nr_skarm+1)), 0.1) - 
  mean(log(red1$nr_skarm+1))

str(D1)

## Define initial values of parameters:
inits <- list(list(alpha.plot = 1,
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

## Old:

## Add distance categories:

# bf$dist_bin[bf$dist <= (50/sqrt(2))] <- "0to35m"
# bf$dist_bin[bf$dist > (50/sqrt(2))] <- "35to50m"
# 
# T1 <- bf
# T1$dist_bin <- "0to50"
# 
# bf <- rbind(bf, T1)



## Add species numbers;
# 
# bf <- as.data.table(bf)
# 
# ## per year/plot/visit;
# bf[, "species_nr" := length(unique(species)),
#    by = c("plot", "visit", "obs_year", "dist_bin")]
# 
# ## and per year/plot/visit/functional group;
# bf[, "species_nr" := length(unique(species)),
#    by = c("plot", "visit", "obs_year", "dist_bin", "func_group")]