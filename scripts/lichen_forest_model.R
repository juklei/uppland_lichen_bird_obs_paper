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
library(ggplot2)
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

## Reduce trees of unknown species:
l_obs <- l_obs[!is.na(l_obs$Tree.species),]

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
l_occ <- melt(l_obs, 
              id.vars = colnames(l_obs)[1:11],
              measure.vars = colnames(l_obs)[c(12:17, 19:100, 102:107)],
              value.name = "observed",
              variable.name = "species")

## Replace NA with 0 in species list to get presence absence per tree:
l_occ$observed[is.na(l_occ$observed)] <- 0

## Exclude combinations of tree species, stem/branches and lichen species that 
## are never observed:
l_occ <- as.data.table(l_occ)
l_occ[, "observable" := ifelse(sum(observed) > 0, "yes", "no"),
      by = c("Tree.species", "Stem.S.Branches.B.T.Total", "species")]

## Export non-observable and let Göran check. Adjust yes/no manually at this 
## point and reimport data set:

# write
# ...
# read

## Exclude non-observable:
l_occ <- l_occ[l_occ$observable == "yes", ]

## Add total observed per tree:

l_temp <- l_occ[, "observed_t" := ifelse(sum(observed) > 0, 1, 0),
                by = c("plot", "Tree.no", "species")]

l_temp$observed <- l_temp$observed_t
l_temp$Stem.S.Branches.B.T.Total <- "T"
l_temp <- unique(l_temp)

l_occ <- rbind(l_occ, l_temp)

l_occ <- l_occ[, -15]

## Bin all deciduous tree species:
levels(l_occ$Tree.species)[levels(l_occ$Tree.species) %in% 
                             c("Qr", "Pt", "Ag", "Bp")] <- "Dc"

## Merge l_occ with forest data:
lof <- merge(l_occ, f_plot[, c(1:3, 6:8, 11:20)], all.x = TRUE, by = "plot") 

## Calculate per tree, which percentage of possible species and the absolute 
## number of species present on that tree:
lof[, c("perc_obs", "richness") := list((sum(observed)/length(observed)),
                                        sum(observed)), 
    by = c("plot", "Tree.no", "Stem.S.Branches.B.T.Total")]

## Reduce to unique rows for percentage occupancy:
lof <- unique(lof[, -12])

## 5. Explore data graphically -------------------------------------------------

g1 <- ggplot(lof, aes(x = PercentBelow5m, 
                      y = perc_obs,
                      fill = Tree.species, 
                      color = Tree.species))
g2 <- geom_point()
g3 <- stat_smooth(method = "lm", size = 2, formula = y ~ log(x))# + I(x^2))            
g4 <- facet_grid(. ~ Stem.S.Branches.B.T.Total) 

g1+g3+g4

## 6. Start analysing lpof ----------------------------------------------------- 

## ...

## -------------------------------END-------------------------------------------
