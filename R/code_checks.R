# rm(list=ls())
# library(readr)
# library(Formula)
# library(dplyr)
# ?distinct
# baltim <- read_csv("/Users/gpiras/SpatialRegime/baltimore/baltim.csv")
#
#
#
# #spatial weights
# crd <- cbind(baltim$X, baltim$Y)
# knn4 <- spdep::knearneigh(crd, k=4)
# wmat <- sphet::listw2dgCMatrix(spdep::nb2listw(spdep::knn2nb(knn4)))
# split = ~ baltim$CITCOU
# form = PRICE  ~ NROOM + NBATH + PATIO + FIREPL +AC + GAR + AGE +LOTSZ +SQFT
# form = PRICE  ~ NROOM + NBATH
# formula <- form
#mod <- regimes(form, baltim, split, vc = "homoskedastic")
