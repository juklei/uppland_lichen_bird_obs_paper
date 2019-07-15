## bird and lichen richness modeled simultaniously with covariates
## 
## First edit: 20190326
## Last edit: 20190715
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
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

bpo <- read.csv("clean/bpo_double_50.csv")

dir("clean")

load("clean/rb_2017.rdata")
bpr_2017 <- zj_pred
load("clean/rb_2018.rdata")
bpr_2018 <- zj_pred
load("clean/lichen_richness.rdata")

## 4. The model ----------------------------------------------------------------

## Create data set with all needed vars:

l_2018 <- data.frame(r_mean_sc = apply(zj_lichen$scaled_rl, 2, mean),
                     r_sd_sc = apply(zj_lichen$scaled_rl, 2, sd),
                     r_mean = apply(zj_lichen$plot_richness, 2, mean),
                     r_sd = apply(zj_lichen$plot_richness, 2, sd),
                     obs_year = 2018,
                     organsim = "lichen",
                     plot = zj_lichen$plotnames)

# b_2017 <- data.frame(r_mean_sc = summary(bpr_2017$scaled_rb, mean)$stat,
#                      r_sd_sc = summary(bpr_2017$scaled_rb, sd)$stat,
#                      r_mean = summary(bpr_2017$richness, mean)$stat,
#                      r_sd = summary(bpr_2017$richness, sd)$stat,
#                      obs_year = 2017,
#                      organsim = "bird",
#                      plot = bpr_2017$plotnames)
# d_all <- rbind(b_2017, l_2018)

b_2018 <- data.frame(r_mean_sc = summary(bpr_2018$scaled_rb, mean)$stat,
                     r_sd_sc = summary(bpr_2018$scaled_rb, sd)$stat,
                     r_mean = summary(bpr_2018$richness, mean)$stat,
                     r_sd = summary(bpr_2018$richness, sd)$stat,
                     obs_year = 2018,
                     organsim = "bird",
                     plot = bpr_2018$plotnames)
l_2018 <- droplevels(l_2018[l_2018$plot %in% b_2018$plot, ])
d_all <- rbind(b_2018, l_2018)

d_all <- merge(d_all, unique(bpo[, c(2,8:13,18)]), all.x = TRUE, by = "plot")

## Export correlations:
d_cor <- bpo[, c(8:13,18)]
d_cor$p_spruce <- bpo$nr_gran/bpo$nr_all_alive
d_cor$p_pine <- bpo$nr_tall/bpo$nr_all_alive
d_cor$p_dec <- bpo$nr_lov/bpo$nr_all_alive
d_cor <- unique(d_cor)
capture.output(cor(d_cor)) %>% write(., "results/correlations.txt")

## Create model data set:
data <- list(nobs = nrow(d_all),
             org = ifelse(d_all$organsim == "lichen", 1, 0),
             r_mean = d_all$r_mean_sc,
             r_sd = d_all$r_sd_sc,
             cd = scale(d_all$PercentAbove3m),
             ud = scale(d_all$PercentBelow3m),
             stand_dbh = scale(d_all$average_dbh_all_alive))

## Add prediction data:

## Understorey density:
data$ud_pred <- seq(min(data$ud), max(data$ud), 0.05)

## Canopy density:
data$cd_pred <- seq(min(data$cd), max(data$cd), 0.05)

## Stand dbh:
data$sdbh_pred <- seq(min(data$stand_dbh), max(data$stand_dbh), 0.05)

str(data)

inits <-  list(list(alpha = 5,
                    beta_org = 0,
                    beta_ud = 0.5,
                    beta_cd = 0.5,
                    beta_sdbh = 0.5,
                    int_ud = 0.1,
                    int_cd = 0.2,
                    int_sdbh = 0.1,
                    plot_sd = 0.1
                    ),
               list(alpha = -2,
                    beta_org = 0.2,
                    beta_ud = 0.7,
                    beta_cd = 0.1,
                    beta_sdbh = 0,
                    int_ud = 0,
                    int_cd = 0,
                    int_sdbh = 0,
                    plot_sd = 0.3
               ),
               list(alpha = 0,
                    beta_org = 0,
                    beta_ud = 1,
                    beta_cd = 1,
                    beta_sdbh = -0.5,
                    int_ud = -0.1,
                    int_cd = -0.2,
                    int_sdbh = -0.1,
                    plot_sd = 0.5
               )
               )

model <- "scripts/JAGS/lichen&bird_JAGS_combined.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 3) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 10000
n.thin <- 5

zc <- coda.samples(jm,
                   variable.names = c("alpha",
                                      "beta_org",
                                      "beta_ud",
                                      "beta_cd",
                                      "beta_sdbh",
                                      "int_ud",
                                      "int_cd",
                                      "int_sdbh",
                                      "plot_sd"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lb_combined_2018_3m_sc.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_lb_combined_2018_3m_sc.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lb_combined_2018_3m_sc.txt")

## 6. Produce and export data for fancy figures --------------------------------

zj_pred <- jags.samples(jm,
                        variable.names = c("rb_ud", 
                                           "rl_ud", 
                                           "rb_cd", 
                                           "rl_cd"
                                           ),
                        n.iter = samples,
                        thin = n.thin)

zj_pred_2017 <- zj_pred
zj_pred_2017$ud <- data$ud_pred#backscale(data$ud_pred, data$ud)
zj_pred_2017$cd <- data$cd_pred#backscale(data$cd_pred, data$cd)
save(zj_pred_2017, file = "clean/combined_pred_2017_7m_sc.rdata")

# zj_pred_sdbh <- zj_pred
# zj_pred_2017_sdbh <- backscale(data$sdbh_pred, data$stand_dbh)
# save(zj_pred_2017_sdbh, file = "clean/combined_pred_2017_sdbh.rdata")

## -------------------------------END-------------------------------------------
