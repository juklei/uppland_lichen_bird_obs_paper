## Make nice figures for the lichen and bird comparison.
## 
## First edit: 20190625
## Last edit: 20200612
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")
require("scales")

## 2. Load and explore data ----------------------------------------------------

## Forest data:
forest <- read.csv("data/forest_data.csv")

## Combined model data:
load("clean/combined_pred_2017_3m_sc.rdata")
c_2017_3 <- zj_pred_2017
load("clean/combined_pred_2017_5m_sc.rdata")
c_2017_5 <- zj_pred_2017
load("clean/combined_pred_2017_7m_sc.rdata")
c_2017_7 <- zj_pred_2017

## Summed model data:
load("clean/summed_pred_2017_3m.rdata")
s_2017_3 <- zj_pred_2017_sum
load("clean/summed_pred_2017_5m.rdata")
s_2017_5 <- zj_pred_2017_sum
load("clean/summed_pred_2017_7m.rdata")
s_2017_7 <- zj_pred_2017_sum

## Taxa model data:
load("clean/rl.rdata")
load("clean/rb_2017.rdata")
bpr_2017 <- zj_bird
load("clean/rb_2018.rdata")
bpr_2018 <- zj_bird

## 3. Make figure plotting scaled against real values to give the readar an 
## impression of magnitude of species richness. Additionally, calculate the 
## formula translating scaled to absolute values:

## For birds:

b_scaled <- summary(bpr_2017$scaled_rb, mean)$stat
b_unscaled <- summary(bpr_2017$richness, mean)$stat

pdf("figures/Figure_S.3a.pdf")
plot(b_scaled, b_unscaled,
     xlab = "Scaled expected bird richness",
     ylab = "Expected bird richness")
dev.off()

## For lichens:

l_scaled <- apply(zj_lichen$scaled_rl, 2, mean)
l_unscaled <- apply(zj_lichen$plot_richness, 2, mean)

pdf("figures/Figure_S.3b.pdf")
plot(l_scaled, l_unscaled,
     xlab = "Scaled expected epiphytic lichen richness",
     ylab = "Expected epiphytic lichen richness")
dev.off()

## Make linear models to extract the correlation function for Fig. 2 between:
## 1) scaled to unscaled for birds:
summary(lm(b_unscaled ~ b_scaled))$coefficients[, "Estimate"]
bs_bird <- function(x) 13.936914 + 2.437788*x 

## 2) scaled to unscaled for lichens
summary(lm(l_unscaled ~ l_scaled))$coefficients[, "Estimate"]
bs_lichen <- function(x) 34.96921 + 6.87766*x

## 3) unscaled birds to unscaled lichens
summary(lm(bs_lichen(0:100) ~ bs_bird(0:100)))$coefficients[, "Estimate"]
bird_to_lichen <- function(x) -4.350598 + 2.821271*x

## 4. Make graphs for birds and lichens for all heightbreaks and both ----------
##    understory and overstory density for 2017

## Make data set:
y_lb <- cbind(summary(c_2017_3$rb_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rb_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rb_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_3$rb_od, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rb_od, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rb_od, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_3$rl_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rl_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rl_ud, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_3$rl_od, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_5$rl_od, quantile, c(.025,.5,.975))$stat,
              summary(c_2017_7$rl_od, quantile, c(.025,.5,.975))$stat)
      
d_lb <- data.frame("richness" = t(y_lb)[,2],
                   "lower" = t(y_lb)[,1],
                   "upper" = t(y_lb)[,3],
                   "density" = c(c_2017_3$ud, c_2017_5$ud, c_2017_7$ud,
                                 c_2017_3$od, c_2017_5$od, c_2017_7$od,
                                 c_2017_3$ud, c_2017_5$ud, c_2017_7$ud,
                                 c_2017_3$od, c_2017_5$od, c_2017_7$od),
                   "organism" = c(rep("birds", length(y_lb[1,])/2),  
                                  rep("epiphytic lichens", length(y_lb[1,])/2)),
                   "story" = c(rep("understory", length(c_2017_3$ud)),
                               rep("understory", length(c_2017_5$ud)),
                               rep("understory", length(c_2017_7$ud)),
                               rep("overstory", length(c_2017_3$od)),
                               rep("overstory", length(c_2017_5$od)),
                               rep("overstory", length(c_2017_7$od)),
                               rep("understory", length(c_2017_3$ud)),
                               rep("understory", length(c_2017_5$ud)),
                               rep("understory", length(c_2017_7$ud)),
                               rep("overstory", length(c_2017_3$od)),
                               rep("overstory", length(c_2017_5$od)),
                               rep("overstory", length(c_2017_7$od))),
                   "hb" = c(rep("heightbreak = 3 meters", length(c_2017_3$ud)),
                            rep("heightbreak = 5 meters", length(c_2017_5$ud)),
                            rep("heightbreak = 7 meters", length(c_2017_7$ud)),
                            rep("heightbreak = 3 meters", length(c_2017_3$od)),
                            rep("heightbreak = 5 meters", length(c_2017_5$od)),
                            rep("heightbreak = 7 meters", length(c_2017_7$od)),
                            rep("heightbreak = 3 meters", length(c_2017_3$ud)),
                            rep("heightbreak = 5 meters", length(c_2017_5$ud)),
                            rep("heightbreak = 7 meters", length(c_2017_7$ud)),
                            rep("heightbreak = 3 meters", length(c_2017_3$od)),
                            rep("heightbreak = 5 meters", length(c_2017_5$od)),
                            rep("heightbreak = 7 meters", length(c_2017_7$od))))

## Transfrom richness to absolute values on bird scale:
d_lb[, 1:3] <- bs_bird(d_lb[, 1:3])

g1 <- ggplot(d_lb, aes(x = density, 
                       y = richness, 
                       fill = organism, 
                       color = organism, 
                       linetype = organism))
g2 <- geom_line(size = 3)
g3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)#, colour = NA)
g4 <- facet_wrap(story ~ hb, scales = "free_x")

png("figures/Figure_2.png", 15000/4, 12000/4, "px", res = 600/4)

g1 + g2 + g3 + g4 +
  geom_hline(yintercept = 13.936914, color = "black", linetype = "dashed")  +
  xlab("vegetation density (% laser returns below/above heightbreak)") + 
  scale_y_continuous(name = "bird richness",
                     sec.axis = sec_axis(trans = bird_to_lichen, 
                                         name = "epiphytic lichen richness")) +
  scale_color_manual(breaks = c("birds", "epiphytic lichens"), 
                     values = c("red", "blue")) + 
  scale_fill_manual(breaks = c("birds", "epiphytic lichens"),
                    values = c("red", "blue")) +
  theme_classic(40) +                  
  theme(legend.position = c(0.17, 0.95), 
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.direction = "horizontal")

dev.off()

## 5. Make graphs for birds and lichens for all heightbreaks and both ----------
##    understory and overstory density for 2017 with summed richness

## Make data set:
y_s <- cbind(summary(s_2017_3$r_ud, quantile, c(.025,.5,.975))$stat,
             summary(s_2017_5$r_ud, quantile, c(.025,.5,.975))$stat,
             summary(s_2017_7$r_ud, quantile, c(.025,.5,.975))$stat,
             summary(s_2017_3$r_od, quantile, c(.025,.5,.975))$stat,
             summary(s_2017_5$r_od, quantile, c(.025,.5,.975))$stat,
             summary(s_2017_7$r_od, quantile, c(.025,.5,.975))$stat)

d_s <- data.frame("richness" = t(y_s)[,2], 
                  "lower" = t(y_s)[,1], 
                  "upper" = t(y_s)[,3],
                  "density" = c(s_2017_3$ud, s_2017_5$ud, s_2017_7$ud,
                                s_2017_3$od, s_2017_5$od, s_2017_7$od),
                  "story" = c(rep("understory", length(s_2017_3$ud)),
                              rep("understory", length(s_2017_5$ud)),
                              rep("understory", length(s_2017_7$ud)),
                              rep("overstory", length(s_2017_3$od)),
                              rep("overstory", length(s_2017_5$od)),
                              rep("overstory", length(s_2017_7$od))),
                  "hb" = c(rep("heightbreak = 3 meters", length(s_2017_3$ud)),
                           rep("heightbreak = 5 meters", length(s_2017_5$ud)),
                           rep("heightbreak = 7 meters", length(s_2017_7$ud)),
                           rep("heightbreak = 3 meters", length(s_2017_3$od)),
                           rep("heightbreak = 5 meters", length(s_2017_5$od)),
                           rep("heightbreak = 7 meters", length(s_2017_7$od))))

p1 <- ggplot(d_s, aes(x = density, y = richness))
p2 <- geom_line(size = 3)
p3 <- geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3)#, colour = NA)
p4 <- facet_wrap(story ~ hb, scales = "free_x")

png("figures/Figure_3.png", 15000/4, 12000/4, "px", res = 600/4)

p1 + p2 + p3 + p4 +
  # geom_vline(xintercept = 0, color = "black", linetype = "dashed")  + 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")  +
  ylab("summed scaled species richness") + 
  xlab("vegetation density (% laser returns below/above heightbreak)") +
  theme_classic(40) 

dev.off()

## 6. Make figure showing Overstory densities vs. the basal area. --------------

## Add basal area: 
forest$ba_tree <- (forest$DBH_ground/2)^2*pi
forest$ba_tree <- forest$ba_tree/10000 ## Calculate m2 from cm2 
forest$ba <- forest$ba_tree*forest$no_all_alive/(pi*.03) ## From subplot to ha

## Create groups for all three heightbreaks and both ud and od:

od_data <- melt(forest[, c(2:4, 16)], measure.vars = colnames(forest[, 2:4]))
levels(od_data$variable) <- c("heightbreak = 3m", 
                              "heightbreak = 5m", 
                              "heightbreak = 7m")
od_data$story <- "overstory"

ud_data <- melt(forest[, c(5:7, 16)], measure.vars = colnames(forest[, 5:7]))
levels(ud_data$variable) <- c("heightbreak = 3m", 
                              "heightbreak = 5m", 
                              "heightbreak = 7m")
ud_data$story <- "understory"

p_data <- rbind(od_data, ud_data)
colnames(p_data)[2] <- "heightbreak"

O <- ggplot(p_data, aes(x = ba, y = value, colour = heightbreak)) +
  geom_rect(aes(xmin = 12, xmax = 20, ymax = Inf, ymin = -Inf, 
                fill = "basal area after thinning"),
            color = NA) +
  scale_fill_manual("", 
                    breaks = "basal area after thinning", 
                    values = "lightgrey") +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              se = FALSE,
              size = 2,
              linetype = "dashed") +
  geom_point(size = 6) +
  facet_grid(story ~ ., scales = "free_y") +
  xlab("basal area (m2/ha)") + 
  ylab("vegetation density (% laser returns above/below heightbreak)") + 
  theme_classic(40) +
  theme(legend.position = c(0.65, 0.5), 
        legend.key.size = unit(3, 'lines'),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "black"))

png("figures/Figure_S.5.png", 9000/4, 12000/4, "px", res = 600/4); O; dev.off()

## 7. Make graphs for difference in species richness per plot ------------------

b_2017 <- data.frame(rb_mean = summary(bpr_2017$richness, mean)$stat,
                     rb_sd = summary(bpr_2017$richness, sd)$stat,
                     plot = bpr_2017$plotnames)
b_2018 <- data.frame(rb_mean = summary(bpr_2018$richness, mean)$stat,
                     rb_sd = summary(bpr_2018$richness, sd)$stat,
                     plot = bpr_2018$plotnames)
b_both <- rbind(b_2017, b_2018)

l_2018 <- data.frame(rl_mean = apply(zj_lichen$plot_richness, 2, mean),
                     rl_sd = apply(zj_lichen$plot_richness, 2, sd),
                     plot = zj_lichen$plotnames)

d_bl <- merge(b_both, l_2018, all.x = TRUE, by = "plot")

## Pearson cc:
cor(d_bl$rl_mean, d_bl$rb_mean)

w1 <- ggplot(d_bl, aes(x = rl_mean, y = rb_mean))
w2 <- geom_point(size = 3)
w4 <- geom_errorbar(aes(ymin = rb_mean - rb_sd, ymax = rb_mean + rb_sd), 
                    color = "red", alpha = 0.5) 
w5 <- geom_errorbarh(aes(xmin = rl_mean - rl_sd, xmax = rl_mean + rl_sd), 
                     color = "blue", alpha = 0.5)
w6 <- annotate("text", 25, 18, label = "pcc = .12", size = 10)

png("figures/Figure_S.4.png", 10000/4, 6500/4, "px", res = 600/4)

w1 + w4 + w5 + w2 + w6 +
  ylab("bird richness") + 
  xlab("epiphytic lichen richness") + 
  theme_classic(40)

dev.off()

## -------------------------------END-------------------------------------------
