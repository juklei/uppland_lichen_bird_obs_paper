## Create data sets to use for Bayesian models in JAGS defined in a different 
## script.
##
## First edit: 20190306
## Last edit: 20190306
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(corrgram)
library(reshape)

## 2. Load and explore data ----------------------------------------------------

dir("data")

occ <- read.csv("data/occ_2016to2018.csv")
f_plot <- read.csv("data/forest_data_uppland_plot.csv")
l_obs <- read.csv("data/lavar_uppland.csv")
exclude <- read.csv("data/non_experiment_include.csv")
ldm <- read.csv("data/long_distance_migrants.csv")

head(occ); head(f_plot); head(l_obs); head(exclude); head(ldm)

## Look at forest explanatory variable correlations and export:

png("figures/forest_correlations_plot.png", 3000, 1000, "px")

corrgram(unique(na.omit(f_plot[, c(6:8, 11:22)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

## 3. Rearrange bird observation data to fit the model presented in the --------
##    AHM book on p. 662. Instead of J indexed on site we want to index it on 
##    the species depending on wether it is long distance migrant or not.

## Exclude all plots in occ which were no lichen inventory occured:
occ <- occ[occ$plot %in% paste0("plot_", unique(l_obs$Plot.no.)), ]

## Exclude plots on which surveys were performed after treatments:
occ <- merge(occ, exclude, by = "plot")
occ <- occ[!(occ$obs_year == 2018 & occ$include_2018 == "N"), ]

## Exclude predators and siskins from obs
occ <- occ[!occ$species %in% c("grona", 
                               "bergk", 
                               "ekore", 
                               "mard", 
                               "gravg"), ]

## Chose relevant columns for this analysis:
occ <- occ[, c("block", "plot", "visit", "obs_year", "obs_time", "species")]

## Make data set so we have all possible combinations of all visits
## and all species seen during the whole survey.
b_occ <- expand.grid.df(unique(occ[, c("block", 
                                       "plot",
                                       "visit",
                                       "obs_year", 
                                       "obs_time")]), 
                        as.data.frame(unique(occ$species)))
colnames(b_occ)[length(b_occ)] <- "species"

## Merge with occ data set to aquire all information of observations
## If a speces was not observed during a visit make obs 0:

## Add obs = 1 to occ:
occ$observed <- 1

## Merge all b_occ with matching b_obs:
b_occ <- merge(b_occ, occ, all.x = TRUE, by = colnames(b_occ))

## NA in observed are actually non-observations, so they will become 0:
b_occ$observed[is.na(b_occ$observed)] <- 0

## If everything went well the number of positive observations should be the 
## nrow of b_obs:
if(sum(b_occ$observed == 1) != nrow(occ)) print("You made a mistake!")

## All long distance migrants could be seen during at ´least three visits. Some
## were seen before, already during the second one. In that case all plots in 
## block in which that species was seen get one more visit for that species.

b_occ$n_visits <- ifelse(b_occ$species %in% ldm$species, 3, 5)

b_occ <- as.data.table(b_occ)




## Has the species been seen in block i during the second visit?
b_occ[b_occ$visit == "second", 
      list("observable_second" = ifelse(sum(observed) > 0, 1, 0)), 
      by = c("obs_year", "block", "species")]

unique(b_occ[b_occ$visit == "second", c("plot", "species")])

## Count during how many visits a species was possible to be seen:
b_occ[b_occ$observable == 1, 
      "n_visits" := nrow(.SD), 
      by = c("obs_year", "plot", "species")]


















## Merge with forest data:
bf_occ <- merge(b_occ, f_plot, all.x = TRUE, by = c("plot", "block"))



## -------------------------------END-------------------------------------------
