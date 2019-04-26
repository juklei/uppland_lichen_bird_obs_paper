## model bird richness in a hierarchical model 
## 
## First edit: 20190409
## Last edit: 20190409
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(coda)
library(magrittr)
library(reshape2)

## 2. Define or source functions used in this script ---------------------------

dir.create("results")
dir.create("results/bpo")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

dir("clean")

bpo <- read.csv("clean/bpo_50.csv")
head(bpo)
str(bpo)

## 4. The model ----------------------------------------------------------------

## Calculate richness and total observation time per plot and year:
bpo <- as.data.table(bpo)
bpo$seen <- ifelse(bpo$n_obs > 0, 1, 0)
bpo[, c("richness", "obs_time_total") := 
      list(sum(seen), sum(obs_time)*max(n_visits)/nlevels(species)), 
    by = c("plot", "obs_year")]

## Create data arrays:

obs <- acast(bpo[, c("plot", "species", "n_obs", "obs_year")],
               formula = plot ~ species ~ obs_year, 
               value.var = "n_obs")
obs[obs[] > 0] <- 1

bpo$obs_time_total <- bpo$n_visits*bpo$obs_time/60
obs_time <- acast(bpo[, c("plot", "species", "obs_time_total", "obs_year")],
                  formula = plot ~ species ~ obs_year, 
                  value.var = "obs_time_total")
obs_time[is.na(obs_time)] <- 0

## Create model data set:
data <- list(nyears = length(unique(bpo$obs_year)),
             nsites = nlevels(bpo$plot),
             nspecies = nlevels(bpo$species),
             obs_time = obs_time,
             obs = obs)

str(data)

## Inits for occ_true
T1 <- data$obs
T1 <- ifelse(T1 > 0, 1, 0)

inits <-  list(list(occ_true = T1,
                    param_obs_time = 0.5,
                    p_occ = 0.5)
               )

model <- "scripts/JAGS/bird_JAGS_bpo_obs_time.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  100000

update(jm, n.iter = burn.in) 

samples <- 20000
n.thin <- 10

zc <- coda.samples(jm,
                   variable.names = c("param_obs_time", 
                                      "p_occ"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/bpo/parameters_bird_obs_time.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bird_obs_time.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/bpo/diagnostics_bird_obs_time.txt")

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("r_year"),
                        n.iter = samples, 
                        thin = n.thin)

## 7. Export data for covariance part ------------------------------------------

## Which plot * year combination does actually exist?
T2 <- unique(bpo[, c("obs_year", "plot")])
T2$exists <- 1
T2 <- acast(T2, formula = plot ~ obs_year, value.var = "exists")
T2[is.na(T2)] <- 0

zj_pred$plots_visited <- T2

save(zj_pred, file = "clean/rb_all_50m_obs_time.rdata")

## -------------------------------END-------------------------------------------
