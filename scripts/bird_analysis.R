## model bird richness in a hierarchical model
## Calculate this script for each year seperately
## 
## First edit: 20190326
## Last edit: 20200612
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)
library(magrittr)
library(reshape2)
library(data.table)

## 2. Define or source functions used in this script ---------------------------

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

dir("clean")

bpo <- read.csv("data/bird_data.csv")
head(bpo)
str(bpo)

dir.create("results")
dir.create("results")
dir.create("figures")

## Calculate number of observation per species:
tmp <- data.table(bpo[, c("species", "n_obs")])
tmp <- tmp[, list("sum" = sum(n_obs)), by = "species"]
capture.output(tmp) %>% write(., "results/bird_occurrences.txt")

## 4. The model ----------------------------------------------------------------

## Do everything for one year each:
bpo <- droplevels(bpo[bpo$obs_year == 2018, ]) ## Change 2017/2018 here !!!!!!!!

## Create data arrays:

nvisits <- acast(bpo[, c("plot", "species", "n_visits")],
                 formula = plot ~ species, 
                 value.var = "n_visits")

nseen <- acast(bpo[, c("plot", "species", "n_obs")],
               formula = plot ~ species, 
               value.var = "n_obs")

## Create model data set:
data <- list(nsites = nlevels(bpo$plot),
             nspecies = nlevels(bpo$species),
             nvisits = nvisits,
             nseen = nseen,
             R = matrix(c(5,0,0,1), ncol = 2),
             df = 3)

str(data)

## Inits for occ_true
T1 <- data$nseen
T1 <- ifelse(T1 > 0, 1, 0)

inits <-  list(list(occ_true = T1,
                    Omega = matrix(c(1,0,0,1), ncol = 2),
                    eta = matrix(0, nrow = data$nspecies, ncol = 2))
               )

model <- "scripts/JAGS/bird_part.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  50000

update(jm, n.iter = burn.in) 

samples <- 100000
n.thin <- 50

zc <- coda.samples(jm,
                   variable.names = c("mu_eta", "probs", "Sigma", 
                                      "p_occ", "p_det"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_bird_part_2018.txt") ## Change 2017/2018 here !!!

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_bird_part_2018.pdf") ## Change 2017/2018 here !!!!!!!!!!!!!!!!
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_bird_part_2018.txt") ## Change 2017/2018 here !!

## 6. Produce results and export data for seperate analysis --------------------

## Produce predictions:
zj_bird <- jags.samples(jm, 
                        variable.names = c("richness", "scaled_rb"),
                        n.iter = samples, 
                        thin = n.thin)

zj_bird$plotnames <- levels(bpo$plot)

save(zj_bird, file = "clean/rb_2018.rdata") ## Change 2017/2018 here !!!!!!!!!!!

## Export expected richness quantiles:
capture.output(summary(zj_bird$richness[], quantile),
               summary(zj_bird$scaled_rb[], quantile)) %>% 
  write(., "results/bird_richness_2018.txt") ## Change 2017/2018 here !!!!!!!!!!

## -------------------------------END-------------------------------------------
