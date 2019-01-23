## In this script l_obsen obs data used in the bird & l_obsen paper is modelled
## with different forest variables. A maximum likelyhood and bayesian approach 
## is used. Different responses are modelled from species to community level.
##
## First edit: 20181207
## Last edit: 20190122
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

corrgram(unique(na.omit(f_subplot[f_subplot$buffer == 10, c(6:8, 11:21)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

png("figures/forest_correlations_plot.png", 3000, 1000, "px")

corrgram(unique(na.omit(f_plot[, c(6:8, 11:22)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

## 4. Rearrange lichen data set to merge with forest data for analysis. -------- 

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
to_GT <- l_occ[, list("observable" = ifelse(sum(observed) > 0, "yes", "no")),
               by = c("Tree.species", "Stem.S.Branches.B.T.Total", "species")]

dir.create("temp")
write.csv(to_GT, "temp/to_GT_observable_lichen.csv", row.names = FALSE)
from_GT <- read.csv("temp/from_GT_observable_lichen.csv")

## Keep all observations and non-observations in l_occ that are observable:
l_occ <- merge(l_occ, 
               from_GT[from_GT$observable == "yes", ],
               by = c("Tree.species", "Stem.S.Branches.B.T.Total", "species"))

## We want to define the block species pool and reduce data set to species 
## observable per block:

## We need "block" from the forest data to define block species pool:
l_occ <- merge(l_occ, f_plot[, c("plot", "block")], by = "plot")

l_occ[, "block_pool" := ifelse(sum(observed) > 0, "yes", "no"),
      by = c("block", "species")]

## Exclude non-available species per block:
l_occ <- l_occ[l_occ$block_pool == "yes", ]

## Add total observed per tree:

l_temp <- l_occ[, "observed_t" := ifelse(sum(observed) > 0, 1, 0),
                by = c("plot", "Tree.no", "species")]

l_temp$observed <- l_temp$observed_t
l_temp$Stem.S.Branches.B.T.Total <- "T"
l_temp <- unique(l_temp)

l_occ <- rbind(l_occ, l_temp)
l_occ <- l_occ[, -17]

## Calculate per tree, which percentage of possible species and the absolute 
## number of species present on that tree:
lo_tree <- unique(l_occ[, list("circle_10m" = circle_10m,
                               "tree_sp" = Tree.species,
                               "tree_dbh" = Tree.diameter.130.cm.above.ground,
                               "perc_obs" = sum(observed)/length(observed), 
                               "richness" = sum(observed)),
                        by = c("plot", "Tree.no", "Stem.S.Branches.B.T.Total")])

## Calculate per plot, which percentage of possible species and the absolute 
## number of species present in that plot:

lo_plot <- l_occ[, list("circle_10m" = unique(circle_10m),
                        "nr_trees" = length(unique(Tree.no)), 
                        "average_tree_dbh" = 
  mean(unique(.SD[, 10:11])$Tree.diameter.130.cm.above.ground, na.rm = TRUE), 
                        "nr_tree_sp" = length(unique(Tree.species)),
                        "richness" = 
  nrow(unique(.SD[.SD$observed == 1, "species"])), 
                        "perc_obs" = 
  nrow(unique(.SD[.SD$observed == 1,"species"]))/length(unique(species))),
                 by = "plot"]

## 5. Merge the different lichen data sets with forest data adjusted for ------- 
##    this question.

## From plot level We want to have the average dbh per plot as a proxy for age
## and the laser measurement:
lof_tree <- merge(lo_tree, 
                  f_plot[, c("plot", "average_dbh_all_alive", "laser_mean")],
                  all.x = TRUE,
                  by = "plot")

lof_plot <- merge(lo_plot, 
                  f_plot[, c("plot", "average_dbh_all_alive", "laser_mean")],
                  all.x = TRUE,
                  by = "plot")

## Then we want to add from the subplot level, the lidar measurments:
lof_tree <- merge(lof_tree,
                  f_subplot[, c(1:3, 6:9)],
                  all.x = TRUE,
                  by = c("plot", "circle_10m"),
                  allow.cartesian = TRUE)

## Bin all deciduous tree species into trivial and complex:
levels(lof_tree$tree_sp)[levels(lof_tree$tree_sp) %in% 
                           c("Ag", "Bp")] <- "Dc_trivial"
levels(lof_tree$tree_sp)[levels(lof_tree$tree_sp) %in% 
                           c("Qr", "Pt")] <- "Dc_complex"
levels(lof_tree$tree_sp)[levels(lof_tree$tree_sp) %in% 
                           c("Dc_complex", "Dc_trivial")] <- "Dc"

lof_plot <- merge(lof_plot,
                  f_subplot[, c(1:3, 6:9)],
                  all.x = TRUE,
                  by = c("plot", "circle_10m"),
                  allow.cartesian = TRUE)

## 6. Explore data graphically -------------------------------------------------

## Categorise plot_dbh for age:

T1 <- mean(lof_tree$average_dbh_all_alive, na.rm = TRUE)

lof_tree$plot_dbh <- ifelse(lof_tree$average_dbh_all_alive > T1,
                            "wide plot dbh", 
                            "narrow plot dbh")
lof_plot$plot_dbh <- ifelse(lof_plot$average_dbh_all_alive > T1,
                            "wide plot dbh", 
                            "narrow plot dbh")

## Factorise nr. of tree species:
lof_plot$nr_tree_sp <- as.factor(lof_plot$nr_tree_sp)

g1 <- ggplot(lof_plot, ## Don't forget to reduce data set for non lidar vars
             aes(x = PercentBelow5m, 
                 y = perc_obs, 
                 fill = nr_tree_sp, 
                 color = nr_tree_sp))
g2 <- geom_point()
g3 <- stat_smooth(method = "lm", size = 2, formula = y ~ log(x))            
g4 <- facet_grid(buffer ~ ., scales = "free") 

dir.create("figures")

png("figures/richness_PB5.png", 2000, 1000, "px")

g1+g3+g4+theme_bw(20)

dev.off()

## 7. Start analysing lof_tree ------------------------------------------------- 

## ...

## -------------------------------END-------------------------------------------






































## Old:

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
lof_pp[, c("nr_trees", "average_dbh", "richness") := 
         list(length(unique(Tree.no)), 
              mean(unique(.SD[, c("Tree.no", 
                                  "Tree.diameter.130.cm.above.ground")])
                   $Tree.diameter.130.cm.above.ground, 
                   na.rm = TRUE),
              length(unique(species))),
       by = c("plot", "Tree.species")]

T_conif <- as.data.table(T_conif)
T_conif[, c("nr_trees", "average_dbh", "richness") := 
          list(length(unique(Tree.no)), 
               mean(unique(.SD[, c("Tree.no", 
                                   "Tree.diameter.130.cm.above.ground")])
                    $Tree.diameter.130.cm.above.ground, 
                    na.rm = TRUE),
               length(unique(species))),
        by = "plot"]

T_all <- as.data.table(T_all)
T_all[, c("nr_trees", "average_dbh", "richness") := 
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