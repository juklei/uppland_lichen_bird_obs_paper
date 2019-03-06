## Create data sets to use for Bayesian models in JAGS defined in a different 
## script.
##
## First edit: 20190306
## Last edit: 20190306
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require("data.table")

## 2. Load and explore data ----------------------------------------------------

dir("data")

occ <- read.csv("data/occ_2016to2018.csv")
f_plot <- read.csv("data/forest_data_uppland_plot.csv")

head(occ); head(f_plot)

## Look at forest explanatory variable correlations and export:

png("figures/forest_correlations_plot.png", 3000, 1000, "px")

corrgram(unique(na.omit(f_plot[, c(6:8, 11:22)])), 
         lower.panel = panel.pie, upper.panel = panel.cor,
         cex.labels = 3)

dev.off()

## 3. Rearrange bird observation data to fit the model presented in the --------
##    AHM book on p. 662. Instead of J indexed on site we want to index it on 
##    the species depending on wether it is long distance migrant or not.




## -------------------------------END-------------------------------------------
