# rm(list=ls())
# library(readr)
# library(Formula)
# library(dplyr)
# ?distinct
# baltim <- read_csv("/Users/gpiras/SpatialRegime/baltimore/baltim.csv")
#
#
# #spatial weights
# crd <- cbind(baltim$X, baltim$Y)
# knn4 <- spdep::knearneigh(crd, k=4)
# wmat <- sphet::listw2dgCMatrix(spdep::nb2listw(spdep::knn2nb(knn4)))
# split = baltim$CITCOU
# form = PRICE  ~ NROOM + NBATH + PATIO + FIREPL +AC + GAR + AGE +LOTSZ +SQFT
# Form <- PRICE ~ NROOM + NBATH + PATIO + FIREPL +AC + GAR + AGE +LOTSZ +SQFT |  NROOM + NBATH + PATIO + FIREPL +AC + GAR
# Form1 <- Formula(Form)
#
# spregime <- function(formula, data, splitv, gw = F){
#
# splitvar <- as.matrix(lm(splitv, data, method="model.frame"))
#  sv <- length(unique(splitvar))
#
#  data <- data[order(splitvar),]
#
#    if("Formula" %in% class(formula)){
#
#     mt <- terms(formula, data = data)
#     mf <- model.frame(formula, data = data)
#
#     y <- model.response(mf)
#     ynm <- all.vars(formula)[1]
#
#     x <- model.matrix(mt, mf)
#
#
#     m <- length(formula)[2]
#     mod <- vector(mode = "list", length = m)
#     modnm <- vector(mode = "list", length = m)
#     lxs <- c()
#     for(i in 1:m){
#       mod[[i]] <- list(model.matrix(formula, data = mf, rhs = i))
#       modnm[[i]] <- attr(mod[[i]][[1]], "dimnames")[[2]]
#       lxs <- c(lxs, length(attr(mod[[i]][[1]], "dimnames")[[2]]))
#     }
#     seq1 <- cumsum(lxs)
#     seq0 <- c(1, cumsum(lxs) +1)[-length(seq1)]
#
#     dt <- matrix(unlist(mod),  nrow = nrow(x), ncol = sum(lxs))
#     colnames(dt) <- unlist(modnm)
#
#
#   dupes <- unlist(lapply(modnm, function(x) lapply(modnm, function(y) intersect(x,y))))
#   xdifnm <- unique(dupes[duplicated(dupes)])
#   xfixnm <- dupes[which(!(dupes %in% dupes[duplicated(dupes)]))]
#   xfix <- matrix(dt[,c(xfixnm)], nrow = nrow(x))
#   colnames(xfix) <- xfixnm
#
#   if(ncol(xfix) ==0) stop("Formula not necessary, use formula instead")
#
#   xdif <- matrix(dt[,c(xdifnm)], nrow = nrow(x))
#   colnames(xdif) <- xdifnm
#   xdifr <- matrix(, nrow = nrow(xfix), ncol=ncol(xdif))
#   xxdif <- matrix(, nrow = nrow(xfix), ncol=ncol(xdif)*sv)
#   zero <- 0
#   sq <- 1 : ncol(xdif)
#   cnames <- vector(mode= "character", length = length(sq))
#   nobss <- vector("numeric", length = length(sv))
#
#    for(i in 1:sv){
#      xdifr <- xdif
#      for(j in 1: ncol(xdifr)) xdifr[splitvar != zero, j]  <- 0
#      nobss[i] <- sum(splitvar == zero, na.rm = TRUE)
#      xxdif[,sq] <- xdifr
#      if(any(xdifnm =="(Intercept)")) xdifnm[which(xdifnm=="(Intercept)")] <- "Intepcept"
#      cnames[sq] <- paste(xdifnm,"_",i,sep="")
#      zero <- zero + 1
#      sq <- sq + ncol(xdif)
#    }
#
#     dataset <- as.data.frame(cbind(y, xfix, xxdif))
#     xnames <- c(xfixnm, cnames)
#
#     colnames(dataset) <- c(ynm, xfixnm, cnames)
#
#     if(any(xnames=="(Intercept)")==1){
#       xnames <- xnames[-which(xnames=="(Intercept)")]
#       form <- as.formula(paste(ynm,"~",  paste(xnames, collapse = " + ")))
#     }
#     else  form <- as.formula(paste(ynm,"~",  paste(xnames, collapse = " + "), paste("-1")))
#
#
#     res <- lm(form , data = dataset)
#
#     #print(summary(res))
#      }
#
#  else{
#
#     dataset <- split.data.frame(data, splitvar)
#     res <- lapply(dataset, lm, formula = formula)
#    # for(i in 1:length(unique(splitvar))) print(summary(res[[i]]))
#  }
#
#  class(res) <- c("spatial_regimes")
#  attr(res, "formula") <- formula
#  attr(res, "nobss") <- nobss
#  attr(res, "splits") <- sv
#  attr(res, "dataset") <- dataset
# return(res)
#
#
#
# }
#
#
#
# vc_mat <-function(obj, type){
#
#
#
# }
#
# summary.spatial_regimes<- function(obj, vc_regimes = c("homo", "hetero"){
#
#   # hetero case
#   if(vc_regimes == "hetero"){
#   ns <- attr(obj,"splits")
#   if("Formula" %in% class(attr(obj,"formula"))){
#
#     nobss <- attr(obj,"nobss")
#     ss <- c(0,cumsum(nobss))
#     lnobss <- length(nobss)
#     uh <- residuals(obj)
#     uhl <- vector(mode = "list", length = lnobss)
#     for(i in 2:(lnobss+1))  uhl[[i-1]] <- uh[ss[i-1]:ss[i]]
#
#
#   }
#   else{
#
#     nobss <- attr(obj,"nobss")
#     uh <- vector(mode = "list", length = lnobss)
#     df <- vector(mode = "list", length = lnobss)
#     vc_mat <- vector(mode = "list", length = lnobss)
#     for(i in 1:ns){
#       uh[[i]] <- residuals(obj[[i]])
#       df[[i]] <- obj[[i]]$df.residual
#       vc_mat[[i]] <- crossprod(uh[[i]]) /  df[[i]]
#     }
#
#     sigma <- rep(unlist(vc_mat), times = nobss)
#     fgls <- lm(formula, data = dataset, weights = 1/sigma)
#   }
#
#   #homo case
#   else{
#
#
#   }
#
#
#
#     cat("Results spatial regime models")
#     summary(obj)
#
#     print()
#   }
#
#
# }
#
# data <- baltim
# splitv <- ~ CITCOU
# hybrid <- spregime(Form1, data = baltim, splitv = ~ CITCOU)
# str(hybrid)
# normal <- spregime(form, data = baltim, splitv  = ~ CITCOU)
# summary(normal)
