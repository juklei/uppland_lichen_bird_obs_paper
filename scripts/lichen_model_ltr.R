## model lichen richness in a hierarchical model with stem diameter and tree
## species at the tree level and vegetation density and stand age at the plot
## level. block is the grouping variable
## 
##
## First edit: 20190125
## Last edit: 20190125
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)

## 2. Define or source functions used in this script ---------------------------

#source()

## 3. Load and explore data ----------------------------------------------------

dir("clean")

ltr <- read.csv("clean/ltr_T_10.csv")
head(ltr)
str(ltr)

## 4. The model ----------------------------------------------------------------

## Plot level explanatory variables need to be reduced to unique rows:
plu <- unique(ltr[, c("plot", 
                      "average_dbh_all_alive", 
                      "PercentAbove5m", 
                      "PercentBelow5m")])
nrow(plu)

data <- list(id = ltr$X,
             block = as.numeric(ltr$block),
             nblock = length(unique(ltr$block)),
             plot = as.numeric(ltr$plot),
             nplot = length(unique(ltr$plot)),
             richness = ltr$richness,
             tree_sp_pine = ifelse(ltr$tree_sp == "Ps", 1, 0),
             tree_sp_spruce = ifelse(ltr$tree_sp == "Pa", 1, 0),
             stem_dbh = scale(ltr$tree_dbh),
             stand_dbh = scale(plu$average_dbh_all_alive),
             canopy_density = scale(plu$PercentAbove5m),
             understory_density = scale(plu$PercentBelow5m))

str(data)

inits <-  list(
  list(beta_pine = 0,
       beta_spruce = 0,
       beta_stem_dbh = 0.5,
       alpha_plot = rep(8, data$nplot),
       sigma_plot = 2,
       alpha_plot_mean = 8,
       beta_stand_dbh = 0.2,
       beta_cdens = -0.2,
       beta_udens = -0.2,
       block_effect = rep(1, data$nblock),
       sigma_block = 2
  ))


model <- "scripts/JAGS/lichen_JAGS_ltr.R"

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
                   variable.names = c("beta_pine",
                                      "beta_spruce",
                                      "beta_stem_dbh",
                                      "alpha_plot",
                                      "alpha_plot_mean",
                                      "beta_stand_dbh",
                                      "beta_cdens",
                                      "beta_udens",
                                      "block_effect",
                                      "sigma_block"), 
                   n.iter = samples, 
                   thin = n.thin)

summary(zc)
plot(zc) 

zj <- jags.samples(jm, 
                   variable.names = c("output"), 
                   n.iter = samples, 
                   thin = n.thin)



#diagnostics
raftery.diag(zc)
heidel.diag(zc)
gelman.diag(zc) #needs at least 2 chains


#useful functions
ecdf(zj$output)(0)
HPDinterval(zc, prob=0.95)
pred < -summary(zj$output, quantile, c(.025,.5,.975))$stat

coda.matrix <- as.matrix(zc[[1]])
head(coda.matrix)


## -------------------------------END-------------------------------------------
