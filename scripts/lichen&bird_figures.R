## Make nice figures for the lichen and bird comparison.
## 
## First edit: 20190221 
## Last edit: 
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")

## 2. Define or source functions used in this script ---------------------------

#...

## 3. Load and explore data ----------------------------------------------------

load("clean/lichen_pred.rdata")
load("clean/bird_pred_veg_age.rdata")
load("clean/bird_pred_trees.rdata")

ltr <- read.csv("clean/ltr_T_10.csv")
bpo <- read.csv("clean/bpo_50.csv")
head(ltr)
str(ltr)

## 4. Make graphs for lichen predictions ---------------------------------------

## Group mean plot for tree species:

## Make data set:
y_ls <- cbind(summary(export_l$dec_mean, quantile, c(.025,.5,.975))$stat,
              summary(export_l$pine_mean, quantile, c(.025,.5,.975))$stat,
              summary(export_l$spruce_mean, quantile, c(.025,.5,.975))$stat)
      
d_ls <- data.frame("richness" = t(y_ls)[,2],
                   "lower" = t(y_ls)[,1],
                   "upper" = t(y_ls)[,3],
                   "species" = c("deciduous", "pine", "spruce"))

g1 <- ggplot(d_ls, aes(x = species, y = richness))
g2 <- geom_point(size = 2)
g3 <- geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3)

png("figures/lichen_tsp.png", 10000, 7000, "px", res = 600)

g1 + g2 + g3 + theme_classic(40) + ylab("lichen richness per tree")

dev.off()

## 5. Make graphs for bird predictions -----------------------------------------

y_btsp <- cbind(summary(export_b_trees$r_dec, 
                        quantile, c(.025,.5,.975))$stat,
                summary(export_b_trees$r_spruce, 
                        quantile, c(.025,.5,.975))$stat,
                summary(export_b_trees$r_pine, 
                        quantile, c(.025,.5,.975))$stat,
                summary(export_b_trees$r_umbr, 
                        quantile, c(.025,.5,.975))$stat)

d_btsp <- data.frame("r" = t(y_btsp)[,2],
                     "lower" = t(y_btsp)[,1],
                     "upper" = t(y_btsp)[,3],
                     "nr_trees" = c(export_b_trees$dec_pred, 
                                    export_b_trees$spruce_pred, 
                                    export_b_trees$pine_pred,
                                    export_b_trees$umbr_pred),
                     "tsp" = c(rep("deciduous", 
                                   length(export_b_trees$dec_pred)),
                               rep("spruce", 
                                   length(export_b_trees$spruce_pred)),
                               rep("pine", 
                                   length(export_b_trees$pine_pred)),
                               rep("umbrella", 
                                   length(export_b_trees$umbr_pred))))

p1 <- ggplot(d_btsp, aes(x = nr_trees, y = r, fill = tsp, color = tsp))
p2 <- geom_line(size = 2)
p3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/bird_tsp.png", 10000, 7000, "px", res = 600)

p1 + p2 + p3 +
  ylab("bird richness") + 
  xlab("number of trees") + 
  theme_classic(40) +                  
  theme(legend.direction = "horizontal", legend.position = c(0.35, 0.95),  
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## 6. Make graphs for both -----------------------------------------------------

## Understory density:

y_ud <- cbind(summary(export_l$r_ud, quantile, c(.025,.5,.975))$stat - 1,
              summary(export_b_veg_age$r_ud, 
                      quantile, c(.025,.5,.975))$stat - 1)

d_ud <- data.frame("r" = t(y_ud)[,2],
                   "lower" = t(y_ud)[,1],
                   "upper" = t(y_ud)[,3],
                   "ud" = c(export_l$ud_pred, export_b_veg_age$ud_pred),
                   "organism" = c(rep("lichens", length(export_l$ud_pred)),
                                  rep("birds", 
                                      length(export_b_veg_age$ud_pred))))

p1 <- ggplot(d_both_ud, aes(x = ud, y = r, fill = organism, color = organism))
p2 <- geom_line(size = 2)
p3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/both_ud.png", 10000, 7000, "px", res = 600)

p1 + p2 + p3 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("understory density") + 
  theme_classic(40) +                  
  theme(legend.position = c(0.9, 0.75), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## Canopy density:

y_cd <- cbind(summary(export_l$r_cd, quantile, c(.025,.5,.975))$stat - 1,
              summary(export_b_veg_age$r_cd,
                      quantile, c(.025,.5,.975))$stat - 1)

d_cd <- data.frame("r" = t(y_cd)[,2],
                   "lower" = t(y_cd)[,1],
                   "upper" = t(y_cd)[,3],
                   "cd" = c(export_l$cd_pred, export_b_veg_age$cd_pred),
                   "organism" = c(rep("lichens", length(export_l$cd_pred)),
                                  rep("birds", 
                                      length(export_b_veg_age$cd_pred))))

q1 <- ggplot(d_both_cd, aes(x = cd, y = r, fill = organism, color = organism))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/both_cd.png", 10000, 7000, "px", res = 600)

q1 + q2 + q3 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("canopy density") + 
  theme_classic(40) +                  
  theme(legend.position = c(0.9, 0.75), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## Stand dbh:

y_sdbh <- cbind(summary(export_l$r_stand_dbh, 
                        quantile, c(.025,.5,.975))$stat - 1,
                summary(export_b_veg_age$r_stand_dbh, 
                        quantile, c(.025,.5,.975))$stat - 1)

d_sdbh <- data.frame("r" = t(y_sdbh)[,2],
                   "lower" = t(y_sdbh)[,1],
                   "upper" = t(y_sdbh)[,3],
                   "sdbh" = c(export_l$stand_dbh_pred, 
                              export_b_veg_age$stand_dbh_pred),
                   "organism" = c(rep("lichens", 
                                      length(export_l$stand_dbh_pred)),
                                  rep("birds", 
                                      length(export_b_veg_age$stand_dbh_pred))))

q1 <- ggplot(d_sdbh, aes(x = sdbh, y = r, fill = organism, color = organism))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/both_sdbh.png", 10000, 7000, "px", res = 600)

q1 + q2 + q3 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("stand dbh (stand age)") + 
  theme_classic(40) +                  
  theme(legend.position = c(0.2, 0.75), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## 6. Make graphs for scale difference in species richness per plot ------------

## Bird richness:
bpo <- as.data.table(bpo)
bpo$seen <- ifelse(bpo$n_obs > 0, 1, 0)
b_2017 <- bpo[bpo$obs_year == 2017, 
              list("richness_birds" = sum(seen)),
              by = "plot"]

## Lichen richness:
combined <- merge(ltr[, c("plot", "richness")], b_2017, by = "plot")
colnames(combined)[2] <- c("richness_lichens")

w1 <- ggplot(combined, aes(x = richness_birds, y = richness_lichens))
w2 <- geom_point(size = 2)
w3 <- stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)

png("figures/both.png", 10000, 7000, "px", res = 600)

w1 + w2 + w3 + 
  ylab("lichen richness per tree") + 
  xlab("bird richness per plot") + 
  theme_classic(40)

dev.off()

## -------------------------------END-------------------------------------------
