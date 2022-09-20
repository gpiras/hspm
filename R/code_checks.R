# # rm(list=ls())
# library(readr)
# library(Formula)
# library(dplyr)
# library(sf)
# library(rgeos)
# library(rgdal)
# ?distinct
#  baltim <- read_csv("/users/gpiras/R dev libraries/SpatialRegime/baltimore/baltim.csv")
#
# #
# #
# #
# # #spatial weights
# # crd <- cbind(baltim$X, baltim$Y)
# # knn4 <- spdep::knearneigh(crd, k=4)
# # wmat <- sphet::listw2dgCMatrix(spdep::nb2listw(spdep::knn2nb(knn4)))
# # split = ~ baltim$CITCOU
# # form = PRICE  ~ NROOM + NBATH + PATIO + FIREPL +AC + GAR + AGE +LOTSZ +SQFT
# form =  PRICE  ~  1 | NROOM + NBATH + PATIO + FIREPL +AC + GAR + AGE +LOTSZ +SQFT
# mod <- regimes(form, baltim, split, vc = "groupwise")
# summary(mod)
#linearHypothesis(mod, c("NBATH_1 - NBATH_2=0"))
#first <- c(0,0,0,0,0,0,0,1,0,0,-1,0,0)
#second <- c(0,0,0,0,0,0,0,0,1,0,0,-1,0)
#third <- c(0,0,0,0,0,0,0,0,0,1,0,0,-1)
#lmat <- rbind(first, second,third)
#wald.test(coef(mod), varb = vcov(mod), L = lmat)
