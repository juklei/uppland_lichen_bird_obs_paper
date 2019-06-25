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

#occ <- read.csv("data/occ_2016to2018.csv")
occ <- read.csv("data/occ_double_2017to2018.csv")
f_plot <- read.csv("data/forest_data_uppland_plot.csv")
l_obs <- read.csv("data/Data_lavar_Almunge15_March_2019.csv")
exclude <- read.csv("data/non_experiment_include.csv")
ldm <- read.csv("data/long_distance_migrants.csv")

head(occ); head(f_plot); head(l_obs); head(exclude); head(ldm)

## Look at forest explanatory variable correlations and export:

dir.create("figures")

png("figures/forest_correlations_plot.png", 3000, 1000, "px")

corrgram(unique(na.omit(f_plot[, c(6:13, 16:27)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

## 3. Rearrange bird observation data to fit the model presented in the --------
##    AHM book on p. 662. Instead of J indexed on site we want to index it on 
##    the species depending on wether it is long distance migrant or not.

# ## Exclude plots in occ which were no lichen inventory occured:
# occ <- occ[occ$plot %in% paste0("plot_", unique(l_obs$Plot.no.)), ]

## Exclude plots which are spruce plantations:
occ <- occ[!occ$plot %in% c("plot_30", "plot_116", "plot_117", "plot_119"), ]

## Exclude nature reserves:
occ <- occ[!occ$plot %in% c("plot_201", "plot_202", "plot_203", "plot_204", 
                            "plot_205", "plot_206", "plot_210", "plot_211",
                            "plot_301", "plot_111", "plot_113", "plot_114",
                            "plot_115", "plot_124"), ]

## Exclude plots on which surveys were performed after treatments:
occ <- merge(occ, exclude, by = "plot")
occ <- occ[!(occ$obs_year == 2018 & occ$include_2018 == "N"), ]

## Exclude predators, birds with large hr and passers from obs:
occ <- occ[!occ$species %in% c("bergk",
                               "duvhk",
                               "ormvk",
                               "ekore",
                               "gravg",
                               "grona", 
                               "korp",
                               "mard",
                               "mindb",
                               "spark"), ]

## Chose relevant columns for this analysis:
occ <- droplevels(occ[, c("block", 
                          "plot",
                          "visit", 
                          "obs_year",
                          "obs_time", 
                          "species")])

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

## Merge all b_occ with matching b_occ:
b_occ <- merge(b_occ, occ, all.x = TRUE, by = colnames(b_occ))

## NA in observed are actually non-observations, so they will become 0:
b_occ$observed[is.na(b_occ$observed)] <- 0

## If everything went well the number of positive observations should be the 
## nrow of b_occ:
if(sum(b_occ$observed == 1) != nrow(occ)) print("You made a mistake!")

## All long distance migrants could be seen during at ´least three visits. Some
## were seen before, already during the second one. In that case all plots in 
## block in which that species was seen get one more visit for that species.

b_occ$n_visits <- ifelse(b_occ$species %in% ldm$species, 3, 5)

b_occ <- as.data.table(b_occ)

## Has the ldm species been seen in block i during the second visit?
T1 <- b_occ[b_occ$visit == "second" & b_occ$species %in% ldm$species, 
            list("observable_second" = ifelse(sum(observed) > 0, 1, 0),
                 "plot" = plot), 
            by = c("obs_year", "block", "species")]

## Reduce to actual observations:
T1 <- unique(T1[T1$observable_second == 1, c("obs_year", "plot", "species")])

## Merge b_occ with T1 to add 1 to rows in T1:

T1$add <- 1
b_occ <- merge(b_occ, T1, all.x = TRUE, by = c("obs_year", "plot", "species"))
b_occ$add[is.na(b_occ$add)] <- 0

b_occ$n_visits <- b_occ$n_visits + b_occ$add

## If everything went well the number of lines with 4 visits should be 5 times
## the nrow of T1:
if(sum(b_occ$n_visits == 4) != 5*nrow(T1)) print("You made a mistake!")

## Count the number of times a species was seen per plot and year:
b_occ[, "n_obs" := sum(observed), by = c("obs_year", "plot", "species")]

## Reduce to needed columns:
b_occ <- unique(b_occ[, c("obs_year", 
                          "block", 
                          "plot", 
                          "species", 
                          "n_obs",
                          "n_visits",
                          "obs_time")])

## Merge with forest data and export:

bf_occ <- merge(b_occ, f_plot[, c(1,8:13,16:17,19:20,25:26)], by = "plot")

dir.create("clean")

write.csv(bf_occ, "clean/bpo_50.csv")

## -------------------------------END-------------------------------------------
