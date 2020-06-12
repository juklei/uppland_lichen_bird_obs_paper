## This data set creates species accumulation data
##
## The species accumulation part:
## We want to shuffle the order of the trees on each plot n times.
## For each shuffle we virtually put the species seen on the first tree
## on all the following trees for that shuffle round and by plot
## For that I create a function which does n shuffelings per plot and apply it 
## then to a list.
##
## First edit: 20190605
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

## 2. Load and explore data ----------------------------------------------------

l_obs <- as.data.table(read.csv("data/epiphytic_lichen_data.csv")) ## Replace with link to dryad!

dir.create("results")
dir.create("clean")
dir.create("figures")

## 3. Create list of matrices per plot that goes into the shuffeling ----------- 

## Change names of plot and circles for consistency:
colnames(l_obs)[2] <- "plot"
l_obs$plot <- as.factor(paste0("plot_", l_obs$plot))

## Exclude plot 119, because not used in bird&lichen comparison, 
## because spruce plantation:
l_obs <- droplevels(l_obs[l_obs$plot != "plot_119", ])

## Reduce l_obs to only uncut trees:
l_obs <- l_obs[!is.na(l_obs$Tree.no), c(2:7, 10:14, 17:137)]

## Reduce to trees of known species:
l_obs <- droplevels(l_obs[!is.na(l_obs$Tree.species) | 
                            l_obs$Tree.species != "check", ])

## Merge Physcia adscendens, P. tenella och P. adscendens/P. tenella to
## P. adscendens/P. tenella"
T1 <- l_obs[, c("Physcia.adscendens", 
                "Physcia.tenella",
                "P..adscendens.tenella")]
l_obs$P..adscendens.tenella <- ifelse(rowSums(T1, na.rm = TRUE) > 0, 1, NA)

## Reduce data set to response columns:
lo_red <- l_obs[, -c("subplot.Center.West.East", "Day.of.inventory", 
                     "Inventoried.by", "Registered.by", 
                     "Cut.stumps.within.the.area..Yes.No.", "Tree.species",
                     "Tree.diameter.130.cm.above.ground", 
                     "X..branches.reaching.ground",
                     "Stem.S.Branches.B.T.Total", 
                     "Bryoria.sp..length.longest.thallus..mm.",
                     "Physcia.adscendens","Physcia.tenella",
                     "Usnea.sp..lenghth.longest.thallus..mm.")]

## Sum species obs to total:
lo_red <- lo_red[, lapply(.SD, sum, na.rm = TRUE), by = c("plot", "Tree.no")]
lo_red[, 3:ncol(lo_red)][lo_red[, 3:ncol(lo_red)] > 1] <- 1 

## 4. Create the data set which is needed to fit the accumultation curve and 
##    produce forest data per plot level with the according order

## Split into a list of data frames per plot:
lo_list <- split(lo_red, by = "plot", keep.by = FALSE)

## Turn the elements into matrices without the tree number:
lo_list <- lapply(lo_list, function(x) as.matrix(x[, -1]))

## Delete all species columns which have never been seen on a plot:
lo_list <- lapply(lo_list, function(x) x[, colSums(x) != 0])

## How many shuffelings?
S <- 100

## Define the function:
sac_create <- function(x) {
  
  out <- matrix(NA, max(l_obs$Tree.no), S)
  for(i in 1:S) {
    ## Select random rows in x:
    D <- x[sample(nrow(x)), ]
    for(j in 1:ncol(D)) {D[min(which(D[, j] == 1)):nrow(D), j] <- 1}
    out[1:nrow(D), i] <- rowSums(D)
  }
  return(out)

}

## Apply the function:
l_sac <- lapply(lo_list, sac_create)

## Store order of the list elements for ordering forest data:
plot_order <- names(l_sac)

## Turn l_sac into an array:
sad <- array(unlist(l_sac), dim = c(max(lengths(l_sac)/S), S, length(l_sac)))

## 5. The model ----------------------------------------------------------------

## Calculate the number of trees by plot:
ntree <- apply(sad, 3, function(x) sum(!is.na(x[,2])))

## Create model data set:
data <- list(nrep = dim(sad)[2],
             nplot = dim(sad)[3],
             ntree = ntree,
             obs = sad)

str(data)

## Prepare inits:
plot_richness <- sad[1,,]
plot_richness[] <- 20
sat_speed <- sad[1,,]
sat_speed[] <- 5

inits <- list(list(plot_richness = plot_richness,
                   lambda_rich = rep(10, data$nrep),
                   lambda_sat = rep(5, data$nrep),
                   sat_speed = sat_speed))

model <- "scripts/JAGS/lichen_part.R"

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
                   variable.names = "plot_richness",
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_lichen_part.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_lichen_part.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_lichen_part.txt")

## Produce validation metrics:
zj_val <- jags.samples(jm,
                       variable.names = "obs_pred",
                       n.iter = 1000,
                       thin = 10)

pred <- summary(zj_val$obs_pred, quantile, c(.025,.5,.975))$stat
x = 0:50

median_mean <- apply(pred[2,,,], c(1,3), mean)
ci_max <- apply(pred[3,,,], c(1,3), max)
ci_min <- apply(pred[1,,,], c(1,3), min)

dev.off()

pdf("figures/sim_vs_obs.pdf")

par(mfrow = c(3, 2))

for(i in 1:data$nplot) {
  
  plot(x, c(0, ci_max[,i]), 
       lty = "dashed", 
       col = "blue", 
       xlab = "tree nr", 
       ylab = "richness", 
       typ = "l")
  lines(x, c(0,median_mean[,i]), col = "blue")
  lines(x, c(0, ci_min[,i]), lty = "dashed", col = "blue")
  
  ## Real data:
  points(rep(which(!is.na(sad[,1,i])), dim(sad)[2]),
         na.omit(as.vector(sad[,,i])))
  
}

dev.off()

## 6. Export lichen richness results and data ----------------------------------

zj_lichen <- jags.samples(jm,
                          variable.names = c("plot_richness", "scaled_rl"),
                          n.iter = samples,
                          thin = n.thin)

zj_lichen$plotnames <- plot_order

save(zj_lichen, file = "clean/rl.rdata")

## Export exected richness quantiles:
capture.output(summary(zj_lichen$plot_richness[], quantile),
               summary(zj_lichen$scaled_rl[], quantile)) %>% 
  write(., "results/lichen_richness.txt")

## -------------------------------END-------------------------------------------
