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

data <- list(id = ltr$X,
             block = as.numeric(ltr$block),
             plot = as.numeric(ltr$plot),
             richness = ltr$richness,
             tree_sp_pine = ifelse(ltr$tree_sp == "Ps", 1, 0),
             tree_sp_spruce = ifelse(ltr$tree_sp == "Pa", 1, 0),
             stem_dbh = ltr$tree_dbh,
             stand_dbh = ltr$average_dbh_all_alive,
             canopy_density = ltr$PercentAbove5m,
             understory_density = ltr$PercentBelow5m)

data

inits <-  list(
  list(
    
  ))


model <- "scripts/JAGS/lichen_ltr_JAGS.R"

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
                   variable.names = c("a", "b", "alpha.mu"), 
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
