## Describe the purpose of the script
## 
##
## First edit: 
## Last edit: 
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

if(!require()){install.packages("")}
require()

## 2. Define or source functions used in this script ---------------------------

source()

## 3. Load and explore data ----------------------------------------------------

dir()

## 4. The model ----------------------------------------------------------------

data <-  list(
  
) 
data

inits <-  list(
  list(
    
  ))


model <- "name_of_JAGS_model_file.R"

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
