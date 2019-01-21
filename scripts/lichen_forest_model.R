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
## are never observed. Export non-observable and let Göran check. Adjust yes/no 
## manually at this point and reimport data set:

l_occ <- as.data.table(l_occ)
# to_GT <- l_occ[, list("observable" = ifelse(sum(observed) > 0, "yes", "no")),
#                by = c("Tree.species", "Stem.S.Branches.B.T.Total", "species")]
# 
# dir.create("temp")
# write.csv(to_GT, "temp/to_GT_observable_lichen.csv", row.names = FALSE)
from_GT <- read.csv("temp/from_GT_observable_lichen.csv")

## Keep all observations and non-observations in l_occ that are observable:
l_occ <- merge(l_occ, 
               from_GT[from_GT$observable == "yes", ],
               by = c("Tree.species", "Stem.S.Branches.B.T.Total", "species"))

## Add total observed per tree:

l_temp <- l_occ[, "observed_t" := ifelse(sum(observed) > 0, 1, 0),
                by = c("plot", "Tree.no", "species")]

l_temp$observed <- l_temp$observed_t
l_temp$Stem.S.Branches.B.T.Total <- "T"
l_temp <- unique(l_temp)

l_occ <- rbind(l_occ, l_temp)

l_occ <- l_occ[, -15]

## Bin all deciduous tree species into trivial and complex:
levels(l_occ$Tree.species)[levels(l_occ$Tree.species) %in%
                             c("Ag", "Bp")] <- "Dc_trivial"
levels(l_occ$Tree.species)[levels(l_occ$Tree.species) %in%
                             c("Qr", "Pt")] <- "Dc_complex"
levels(l_occ$Tree.species)[levels(l_occ$Tree.species) %in%
                             c("Dc_complex", "Dc_trivial")] <- "Dc"

## Merge l_occ with forest data:
lof <- merge(l_occ, f_plot[, c(1:3, 6:8, 11:20)], all.x = TRUE, by = "plot") 
# lof <- merge(l_occ, 
#              f_subplot[, c(1:3, 6:8, 11)], 
#              all.x = TRUE, 
#              by = c("plot", "circle_10m")) 

## Calculate per tree, which percentage of possible species and the absolute 
## number of species present on that tree:
lof[, c("perc_obs", "richness") := list((sum(observed)/length(observed)),
                                        sum(observed)), 
    by = c("plot", "Tree.no", "Stem.S.Branches.B.T.Total")]

## Reduce to unique rows, excluding species to look at per tree richness and 
## perc_obs, representing percentage of diversity potential:
lof <- as.data.frame(lof)
lof_pt <- unique(lof[, -which(colnames(lof) == "species")])

## Calculate cumultative lichen diversity per plot for different combinations of
## tree species to look at relative contribution od tree species:

## Only for total per tree observed:
lof_pp <- lof[lof$Stem.S.Branches.B.T.Total == "T" & lof$observed == 1, 
              -which(colnames(lof) %in% c("perc_obs", "richness"))] 

## Add conifer and all as tree species names:

T_conif <- lof_pp[lof$Tree.species %in% c("Pa", "Ps"), ]
T_conif$Tree.species <- "conifer"

T_all <- lof_pp
T_all$Tree.species <- "all"

## Calculate richness, nr. of trees and mean dbh per plot for each tree group:

lof_pp <- as.data.table(lof_pp)
lof_pp[, c("nr_tress", "average_dbh", "richness") := 
         list(length(unique(Tree.no)), 
              mean(unique(.SD[, c("Tree.no", 
                                  "Tree.diameter.130.cm.above.ground")])
                   $Tree.diameter.130.cm.above.ground, 
                   na.rm = TRUE),
              length(unique(species))),
       by = c("plot", "Tree.species")]

T_conif <- as.data.table(T_conif)
T_conif[, c("nr_tress", "average_dbh", "richness") := 
          list(length(unique(Tree.no)), 
               mean(unique(.SD[, c("Tree.no", 
                                   "Tree.diameter.130.cm.above.ground")])
                    $Tree.diameter.130.cm.above.ground, 
                    na.rm = TRUE),
               length(unique(species))),
        by = "plot"]

T_all <- as.data.table(T_all)
T_all[, c("nr_tress", "average_dbh", "richness") := 
        list(length(unique(Tree.no)), 
             mean(unique(.SD[, c("Tree.no", 
                                 "Tree.diameter.130.cm.above.ground")])
                  $Tree.diameter.130.cm.above.ground, 
                  na.rm = TRUE),
             length(unique(species))),
      by = "plot"]


lof_pp <- rbind(lof_pp, T_conif, T_all)

## Chose desired columns:
lof_pp <- unique(lof_pp[, c(1:2, 15, 17:32)])

## 5. Explore data graphically -------------------------------------------------

## Categorise dbh:
T1 <- mean(lof_pt$Tree.diameter.130.cm.above.ground, na.rm = TRUE)
lof_pt$dbh_cat <- ifelse(lof_pt$Tree.diameter.130.cm.above.ground > T1, 
                         "wide stem", 
                         "narrow stem")

## Categorise forest_height:
T2 <- mean(lof_pp$ElevP95, na.rm = TRUE)
lof_pp$E95_cat <- ifelse(lof_pp$ElevP95 > T2, "high forest", "low forest")

g1 <- ggplot(lof_pp, aes(x = PercentBelow5m, 
                         y = richness,
                         fill = Tree.species, 
                         color = Tree.species))
g2 <- geom_point()
g3 <- stat_smooth(method = "lm", size = 2, formula = y ~ log(x))            
g4 <- facet_grid(E95_cat ~ .) 

dir.create("figures")

png("figures/richness_PB.png", 2000, 1000, "px")

g1+g3+g4+theme_bw(30)

dev.off()

## 6. Start analysing lof_pt ----------------------------------------------------- 

## ...

## -------------------------------END-------------------------------------------
