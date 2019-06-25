## model bird richness in a hierarchical model 
## 
## First edit: 20190307
## Last edit: 20190617
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

dir.create("results")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## 3. Load and explore data ----------------------------------------------------

dir("clean")

bpo <- read.csv("clean/bpo_50.csv")
head(bpo)
str(bpo)

## 4. The model ----------------------------------------------------------------

## Calculate richness and total observation time per plot and year:
bpo <- as.data.table(bpo)
bpo$seen <- ifelse(bpo$n_obs > 0, 1, 0) ## At least two obs per plot and season
bpo[, c("richness", "obs_time_total") := 
      list(sum(seen), sum(obs_time)*max(n_visits)/nlevels(species)), 
    by = c("plot", "obs_year")]

## Reduce to unique rows:
bpr <- unique(bpo[, c(2:4,9:20,22:23)])

## Create model data set:
data <- list(nobs = nrow(bpr),
             nsites = nlevels(bpr$plot),
             plot = as.numeric(bpr$plot),
             obs_year = as.numeric(bpr$obs_year),
             years = unique(as.numeric(bpr$obs_year)),
             richness = bpr$richness,
             obs_time = bpo$obs_time)

str(data)

inits <-  list(list(richness_true = data$richness,
                    plot_richness = rep(5, data$nsites),
                    param_obs_time = 5,
                    sigma_year = 2))

model <- "scripts/JAGS/bird_JAGS_bpr_part.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 1000
n.thin <- 5

zc <- coda.samples(jm,
                   variable.names = c("param_obs_time",
                                      "plot_richness",
                                      "sigma_year"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_bpr_part.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bpr_part.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_bpr_part.txt")

# ## Produce validation metrics:
# zj_val <- jags.samples(jm,
#                        variable.names = c("mean_richness",
#                                           "mean_richness_sim",
#                                           "p_mean",
#                                           "cv_richness",
#                                           "cv_richness_sim",
#                                           "p_cv",
#                                           "fit",
#                                           "fit_sim",
#                                           "p_fit"),
#                        n.iter = samples,
#                        thin = n.thin)
# 
# ## Fit of mean:
# plot(zj_val$mean_richness,
#      zj_val$mean_richness_sim,
#      xlab = "mean real",
#      ylab = "mean simulated",
#      cex = .05)
# abline(0, 1)
# p <- summary(zj_val$p_mean, mean)
# text(x = 0.5, y = 0.66, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Fit of variance:
# plot(zj_val$cv_richness,
#      zj_val$cv_richness_sim,
#      xlab = "cv real",
#      ylab = "cv simulated",
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_cv, mean)
# text(x = 1.5, y = 1.9, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Overall fit:
# plot(zj_val$fit,
#      zj_val$fit_sim,
#      xlab = "ssq real",
#      ylab = "ssq simulated",
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_fit, mean)
# text(x = 1300, y = 1200, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Produce and export figures -----------------------------------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("r_ud",
                                           "r_cd",
                                           "r_stand_dbh"),
                        n.iter = samples, 
                        thin = n.thin)

## 7. Export data for fancy figures --------------------------------------------

export_b <- zj_pred
export_b$ud_pred <- data$ud_pred
export_b$cd_pred <- data$cd_pred
export_b$stand_dbh_pred <- data$stand_dbh_pred

save(export_b, file = "clean/bird_pred.rdata")

## -------------------------------END-------------------------------------------
