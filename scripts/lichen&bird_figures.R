## Make nice figures for the lichen and bird comparison.
## 
## First edit: 20190429
## Last edit: 20190507
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
load("clean/bird_pred.rdata")

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

png("figures/lichen_tsp.png", 10000/4, 7000/4, "px", res = 600/4)

g1 + g2 + g3 + theme_classic(40) + ylab("lichen richness per tree")

dev.off()

## 5. Make graphs for both -----------------------------------------------------

## Understory density:

y_ud <- cbind(summary(export_l$r_ud, quantile, c(.025,.5,.975))$stat - 1,
              summary(export_b$r_ud, quantile, c(.025,.5,.975))$stat - 1)

d_ud <- data.frame("r" = t(y_ud)[,2],
                   "lower" = t(y_ud)[,1],
                   "upper" = t(y_ud)[,3],
                   "ud" = c(export_l$ud_pred, export_b$ud_pred),
                   "organism" = c(rep("lichens", length(export_l$ud_pred)),
                                  rep("birds", length(export_b$ud_pred))))

p1 <- ggplot(d_ud, aes(x = ud, y = r, fill = organism, color = organism))
p2 <- geom_line(size = 2)
p3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/both_ud.png", 10000/4, 7000/4, "px", res = 600/4)

p1 + p2 + p3 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("understory density") + 
  scale_color_manual(breaks = c("birds", "lichens"), 
                     values = c("red", "blue")) + 
  scale_fill_manual(breaks = c("birds", "lichens"),
                    values = c("red", "blue")) +
  theme_classic(40) +                  
  theme(legend.position = c(0.9, 0.75), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## Canopy density:

y_cd <- cbind(summary(export_l$r_cd, quantile, c(.025,.5,.975))$stat - 1,
              summary(export_b$r_cd, quantile, c(.025,.5,.975))$stat - 1)

d_cd <- data.frame("r" = t(y_cd)[,2],
                   "lower" = t(y_cd)[,1],
                   "upper" = t(y_cd)[,3],
                   "cd" = c(export_l$cd_pred, export_b$cd_pred),
                   "organism" = c(rep("lichens", length(export_l$cd_pred)),
                                  rep("birds", length(export_b$cd_pred))))

q1 <- ggplot(d_cd, aes(x = cd, y = r, fill = organism, color = organism))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/both_cd.png", 10000/4, 7000/4, "px", res = 600/4)

q1 + q2 + q3 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("canopy density") + 
  scale_color_manual(breaks = c("birds", "lichens"), 
                     values = c("red", "blue")) + 
  scale_fill_manual(breaks = c("birds", "lichens"),
                    values = c("red", "blue")) +
  theme_classic(40) +                  
  theme(legend.position = c(0.9, 0.75), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## ud and cd combined:

d_ud$cat <- "understory - below 5m"
colnames(d_ud)[4] <- "veg_dens"
d_cd$cat <- "canopy - above 5m"
colnames(d_cd)[4] <- "veg_dens"
d_both <- rbind(d_ud, d_cd)

v1 <- ggplot(d_both, aes(x = veg_dens, 
                         y = r, 
                         fill = organism, 
                         color = organism,
                         linetype = cat))
v2 <- geom_line(size = 2)
v3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/both_both.png", 10000/4, 7000/4, "px", res = 600/4)

v1 + v2 + v3 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("vegetation density") + 
  scale_color_manual(breaks = c("birds", "lichens"), 
                     values = c("red", "blue")) + 
  scale_fill_manual(breaks = c("birds", "lichens"),
                    values = c("red", "blue")) +
  theme_classic(40) +                  
  theme(legend.position = c(0.32, 0.9), 
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## Stand dbh:

y_sdbh <- cbind(summary(export_l$r_stand_dbh, 
                        quantile, c(.025,.5,.975))$stat - 1,
                summary(export_b$r_stand_dbh, 
                        quantile, c(.025,.5,.975))$stat - 1)

d_sdbh <- data.frame("r" = t(y_sdbh)[,2],
                   "lower" = t(y_sdbh)[,1],
                   "upper" = t(y_sdbh)[,3],
                   "sdbh" = c(export_l$stand_dbh_pred, export_b$stand_dbh_pred),
                   "organism" = c(rep("lichens", 
                                      length(export_l$stand_dbh_pred)),
                                  rep("birds", 
                                      length(export_b$stand_dbh_pred))))

q1 <- ggplot(d_sdbh, aes(x = sdbh, y = r, fill = organism, color = organism))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)

png("figures/both_sdbh.png", 10000/4, 7000/4, "px", res = 600/4)

q1 + q2 + q3 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("stand dbh (stand age)")  + 
  scale_color_manual(breaks = c("birds", "lichens"), 
                     values = c("red", "blue")) + 
  scale_fill_manual(breaks = c("birds", "lichens"),
                    values = c("red", "blue")) +
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

## Pearson cc:
combined <- as.data.table(combined)
pcc <- cor(combined[, list("mean_richness_lichens" = mean(richness_lichens), 
                           "richness_birds" = unique(richness_birds)), 
                    by = plot][, 2:3])

w1 <- ggplot(combined, aes(x = richness_birds, y = richness_lichens))
w2 <- geom_point(size = 2)
w3 <- stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)
w4 <- annotate("text", 8, 20, label = "pcc = -0.012", size = 10)

png("figures/both.png", 10000/4, 7000/4, "px", res = 600/4)

w1 + w2 + w4 +
  ylab("lichen richness per tree") + 
  xlab("bird richness per plot") + 
  theme_classic(40)

dev.off()

## -------------------------------END-------------------------------------------
