## In this script l_obsen obs data used in the bird & l_obsen paper is modelled
## with different forest variables. A maximum likelyhood and bayesian approach 
## is used. Different responses are modelled from species to community level.
##
## First edit: 20181207
## Last edit: 20190116
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(dplyr)
library(boot)
library(corrgram)
# library(rjags)
# library(coda)

## 2. Define or source functions used in this script ---------------------------

#source()

## 3. Load and explore data ----------------------------------------------------

dir("data")
f_subplot <- read.csv("data/forest_data_uppland_subplot.csv")
f_plot <- read.csv("data/forest_data_uppland_plot.csv")
l_obs <- read.csv("data/lavar_uppland.csv")

head(f_subplot); head(f_plot); head(l_obs)

## Look at forest explanatory variable correlations and export:

dir.create("figures")

png("figures/forest_correlations_subplot.png", 3000, 1000, "px")

corrgram(unique(na.omit(f_subplot[, c(6:8, 11)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

png("figures/forest_correlations_plot.png", 3000, 1000, "px")

corrgram(unique(na.omit(f_plot[, c(6:8, 11:20)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

## 4. Rearrange lichen data set and merge with forest data for analysis. ------- 

## Change names of plot and circles for merging with forest data:

colnames(l_obs)[2] <- "plot"
l_obs$plot <- as.factor(paste0("plot_", l_obs$plot))

colnames(l_obs)[3] <- "circle_10m"
levels(l_obs$circle_10m) <- c("middle", "east", "west")

## Exclude cut stumps and make seperate data set:
l_obs_cut <- l_obs[is.na(l_obs$Tree.no), c(2:3, 8:9)]

## Reduce l_obs to only uncut trees:
l_obs <- l_obs[!is.na(l_obs$Tree.no), c(2:7, 10:14, 17:112)]

## Exclude thallus length measurments and make seperate data set:

l_obs_thallus <- melt(l_obs,
                      id.vars = colnames(l_obs)[1:11],
                      measure.vars = colnames(l_obs)[c(18, 101)],
                      value.name = "thallus_length",
                      variable.name = "species")

## Replace NA with 0 in species list to get thallus length 0:
l_obs_thallus$thallus_length[is.na(l_obs_thallus$thallus_length)] <- 0

## Merge l_obs_thallus with forest data:
ltf <- merge(l_obs_thallus, 
             f_subplot[, c(1:3, 6:8, 11)], 
             all.x = TRUE, 
             by = c("plot", "circle_10m")) 

## Make data frame based on occupancy of all species per tree:
l_obs_occ <- melt(l_obs, 
                  id.vars = colnames(l_obs)[1:11],
                  measure.vars = colnames(l_obs)[c(12:17, 19:100, 102:107)],
                  value.name = "observed",
                  variable.name = "species")

## Replace NA with 0 in species list to get presence absence per tree:
l_obs_occ$observed[is.na(l_obs_occ$observed)] <- 0

## Exclude combinations of tree species, stem/branches and lichen species that 
## are never observed:
l_obs_occ <- as.data.table(l_obs_occ)
l_obs_occ[, "observable" := ifelse(sum(observed) > 0, "yes", "no"),
          by = c("Tree.species", "Stem.S.Branches.B.T.Total", "species")]

## Export non-observable and let Göran check. Adjust yes/no manually at this 
## point and reimport data set:

# write
# ...
# read

## Exclude non-observable:
l_obs_occ_red <- l_obs_occ[l_obs_occ$observable == "yes", ]

## Merge l_obs_occ_red with forest data:
lof <- merge(l_obs_occ_red, f_plot[, c(1:3, 6:8, 11:20)], 
             all.x = TRUE, 
             by = "plot") 

## Calculate per tree, which percentage of possible species is present on tree:
lof[, "perc_obs" := (sum(observed)/length(observed)), 
    by = c("plot", "Tree.no", "Stem.S.Branches.B.T.Total")]

## Reduce to unique rows for percentage occupancy:
lpof <- unique(lof[, -12])

## 5. Explore data graphically -------------------------------------------------




## 6. Start analysing lpof ----------------------------------------------------- 

## ...

## -------------------------------END-------------------------------------------
