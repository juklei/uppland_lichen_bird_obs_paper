## Make nice figures for the lichen and bird comparison.
## 
## First edit: 20190625
## Last edit: 20190625
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

load("clean/combined_pred_2017_3m.rdata")
c_2017_3 <- zj_pred_2017
load("clean/combined_pred_2017_5m.rdata")
c_2017_5 <- zj_pred_2017
load("clean/combined_pred_2017_7m.rdata")
c_2017_7 <- zj_pred_2017

load("clean/combined_pred_2017_sdbh.rdata")
str(zj_pred_2017_sdbh)

load("clean/lichen_richness.rdata")
str(zj_lichen)

load("clean/rb_2017.rdata")
bpr_2017 <- zj_pred
str(bpr_2017)

## 4. Make graphs for birds and lichens for all heightbreaks and both ----------
##    understory and canopy density for 2017

## Make data set:
y_lb <- cbind(summary(c_2017_3$rb_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rb_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rb_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_3$rb_cd, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rb_cd, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rb_cd, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_3$rl_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rl_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rl_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_3$rl_cd, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rl_cd, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rl_cd, quantile, c(.025,.5,.975))$stat)
      
d_lb <- data.frame("richness" = t(y_lb)[,2],
                   "lower" = t(y_lb)[,1],
                   "upper" = t(y_lb)[,3],
                   "density" = c(c_2017_3$ud, c_2017_5$ud, c_2017_7$ud,
                                 c_2017_3$cd, c_2017_5$cd, c_2017_7$cd,
                                 c_2017_3$ud, c_2017_5$ud, c_2017_7$ud,
                                 c_2017_3$cd, c_2017_5$cd, c_2017_7$cd),
                   "organism" = c(rep("birds", length(y_lb[1,])/2),  
                                  rep("lichens", length(y_lb[1,])/2)),
                   "story" = c(rep("understory", length(c_2017_3$ud)),
                               rep("understory", length(c_2017_5$ud)),
                               rep("understory", length(c_2017_7$ud)),
                               rep("canopy", length(c_2017_3$cd)),
                               rep("canopy", length(c_2017_5$cd)),
                               rep("canopy", length(c_2017_7$cd)),
                               rep("understory", length(c_2017_3$ud)),
                               rep("understory", length(c_2017_5$ud)),
                               rep("understory", length(c_2017_7$ud)),
                               rep("canopy", length(c_2017_3$cd)),
                               rep("canopy", length(c_2017_5$cd)),
                               rep("canopy", length(c_2017_7$cd))),
                   "hb" = c(rep("heightbreak = 3 meters", length(c_2017_3$ud)),
                            rep("heightbreak = 5 meters", length(c_2017_5$ud)),
                            rep("heightbreak = 7 meters", length(c_2017_7$ud)),
                            rep("heightbreak = 3 meters", length(c_2017_3$cd)),
                            rep("heightbreak = 5 meters", length(c_2017_5$cd)),
                            rep("heightbreak = 7 meters", length(c_2017_7$cd)),
                            rep("heightbreak = 3 meters", length(c_2017_3$ud)),
                            rep("heightbreak = 5 meters", length(c_2017_5$ud)),
                            rep("heightbreak = 7 meters", length(c_2017_7$ud)),
                            rep("heightbreak = 3 meters", length(c_2017_3$cd)),
                            rep("heightbreak = 5 meters", length(c_2017_5$cd)),
                            rep("heightbreak = 7 meters", length(c_2017_7$cd)))
                   )

g1 <- ggplot(d_lb, aes(x = density, 
                       y = richness, 
                       fill = organism, 
                       color = organism, 
                       linetype = story))
g2 <- geom_line(size = 2)
g3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA)
g4 <- facet_grid(. ~ hb)

png("figures/both_both.png", 20000/4, 8000/4, "px", res = 600/4)

g1 + g2 + g3 + g4 +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  + 
  ylab("richness relative to the intercept") + 
  xlab("vegetation density") + 
  scale_color_manual(breaks = c("birds", "lichens"), 
                     values = c("red", "blue")) + 
  scale_fill_manual(breaks = c("birds", "lichens"),
                    values = c("red", "blue")) +
  theme_classic(40) +                  
  theme(legend.position = c(0.35, 0.9), 
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## 5. Make graphs for age with absolute richness values ------------------------

## Stand dbh:

y_sdbh <- cbind(summary(zj_pred_2017_sdbh$rb_sdbh, 
                        quantile, c(.025,.5,.975))$stat,
                summary(zj_pred_2017_sdbh$rl_sdbh, 
                        quantile, c(.025,.5,.975))$stat)

d_sdbh <- data.frame("r" = t(y_sdbh)[,2],
                     "lower" = t(y_sdbh)[,1],
                     "upper" = t(y_sdbh)[,3],
                     "sdbh" = rep(zj_pred_2017_sdbh$sdbh, 2),
                     "organism" = c(rep("birds", length(y_sdbh[1,])/2),
                                    rep("lichens", length(y_sdbh[1,])/2))
                     )

q1 <- ggplot(d_sdbh, aes(x = sdbh, y = r, fill = organism, color = organism))
q2 <- geom_line(size = 2)
q3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, colour = NA)

png("figures/both_sdbh.png", 10000/4, 8000/4, "px", res = 600/4)

q1 + q2 + q3 +
  geom_hline(yintercept = min(d_sdbh$upper[d_sdbh$organism == "lichens"]),
             color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = min(d_sdbh$upper[d_sdbh$organism == "birds"])
             , color = "black", linetype = "dashed")  + 
  ylab("richness") + 
  xlab("stand age (stand dbh in cm)")  + 
  scale_color_manual(breaks = c("birds", "lichens"), 
                     values = c("red", "blue")) + 
  scale_fill_manual(breaks = c("birds", "lichens"),
                    values = c("red", "blue")) +
  theme_classic(40) +                  
  theme(legend.position = c(0.12, 0.8), 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'))

dev.off()

## 6. Make graphs for scale difference in species richness per plot ------------

b_2017 <- data.frame(rb_mean = summary(bpr_2017$richness, mean)$stat,
                     rb_sd = summary(bpr_2017$richness, sd)$stat,
                     plot = bpr_2017$plotnames)

l_2018 <- data.frame(rl_mean = apply(zj_lichen$plot_richness, 2, mean),
                     rl_sd = apply(zj_lichen$plot_richness, 2, sd),
                     plot = zj_lichen$plotnames)

d_bl <- merge(b_2017, l_2018, by = "plot")

## Pearson cc:
pcc <- cor(d_bl$rl_mean, d_bl$rb_mean)

w1 <- ggplot(d_bl, aes(x = rl_mean, y = rb_mean))
w2 <- geom_point(size = 3)
# w3 <- stat_smooth(method = "lm", formula = y ~ poly(x, 2), #se = FALSE, 
#                   color = "black")
w4 <- geom_errorbar(aes(ymin = rb_mean - rb_sd, ymax = rb_mean + rb_sd), 
                    color = "red", alpha = 0.5) 
w5 <- geom_errorbarh(aes(xmin = rl_mean - rl_sd, xmax = rl_mean + rl_sd), 
                     color = "blue", alpha = 0.5)
w6 <- annotate("text", 25, 18, label = "pcc = .096", size = 10)

png("figures/both_cor.png", 10000/4, 6500/4, "px", res = 600/4)

w1 + w2 + w4 + w5 + w6 +
  ylab("bird richness") + 
  xlab("lichen richness") + 
  theme_classic(40)

dev.off()

## -------------------------------END-------------------------------------------
