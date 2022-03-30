##### Functions for spatial lag regimes model####
#' Estimation of spatial regime models
#' @name lag_regimes
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{data.frame}.
#' @param listw a spatial weighting matrix
#' @param rgv variable to identify the regimes
#' @param wy_rg default \code{wy_rg = TRUE}, the lagged dependent variable varies for regime
#' @param het heteroskedastic variance-covariance matrix
#' @param cl record calls
#' @param object an object of class lag_regime
#' @param ... additional arguments
#' @param x an object of class lag_regime
#' @param digits number of digits

#'
lag_regimes <- function(formula, data, listw, rgv,
                   het, cl, wy_rg){




  intro <- iv.lag.data.prep.regimes(formula, data = data, listw = listw,
                           listw2 = NULL, wy_rg = wy_rg, rgv = rgv)


  y        <- intro[[1]]
  Hmat     <- intro[[2]]
  Zmat     <- intro[[3]]
  colnames.end<- intro[[4]]
  colnames.instr <- intro[[5]]



  res <- spatial.ivreg.regimes(as.matrix(y), as.matrix(Zmat), as.matrix(Hmat), het)
  res <- list(res, cl, colnames.end,  colnames.instr)
  class(res) <- "lag_regimes"
  return(res)

}

### S3 methods ----

#' @rdname lag_regimes
#' @method coef lag_regimes
#' @export
coef.lag_regimes <- function(object, ...){
  object[[1]][[1]]
}


#' @rdname lag_regimes
#' @method vcov lag_regimes
#' @import stats
#' @export
vcov.lag_regimes <- function(object, ...){
  V <- object[[1]][[2]]
  return(V)
}


#' @rdname lag_regimes
#' @method print lag_regimes
#' @import stats
#' @export
print.lag_regimes <- function(x,
                            digits = max(3, getOption("digits") - 3),
                            ...)
{
  cat("Call:\n")
  print(x[[2]])

  cat("\nCoefficients:\n")

  print.default(format(drop(coef(x)), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}



#' @rdname lag_regimes
#' @method summary lag_regimes
#' @import stats
#' @export
summary.lag_regimes <- function(object, ...){
  b                   <- coef(object)
  std.err             <- sqrt(diag(vcov(object)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.lag_regimes", "lag_regimes")
  return(object)
}


#' @rdname lag_regimes
#' @method print summary.lag_regimes
#' @import stats
#' @export
print.summary.lag_regimes <- function(x,
                                    digits = max(5, getOption("digits") - 3),
                                    ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                          Spatial Lag Regimes Model \n")
  cat("        ------------------------------------------------------------\n")
  cat("\nCall:\n")
  cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

  cat("\nEndogenous variables:\n")

  cat(paste(unlist(x[[3]]), sep=" "))

  cat("\nInstruments:\n")

  cat(paste(x[[4]], sep=" "))
  invisible(x)
}





iv.lag.data.prep.regimes <- function(formula, data, rgv, listw, listw2, wy_rg ){


  #define the split variable
  if(is.null(rgv)) stop("regimes variable not specified")
  if(class(rgv) != "formula") stop("regimes variable has to be a formula")

  n                <- dim(data)[1]
  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- as.numeric(unique(splitvar))
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))

#define w matrix
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
  if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)
  if(inherits(listw,"matrix"))  Ws <- Matrix(listw)
  if(inherits(listw,"Matrix"))  Ws <- listw

  if (n != nrow(Ws))
    stop("Input data and weights have different dimension")

 ####### extract the elements of the formula
  F1 <- Formula(formula)
  parts <- length(F1)

  if(parts[2] != 6) stop("Formula should have six parts")

  mf <- model.frame(F1, data = data)
  y <- model.part(F1, data = mf, lhs = 1, drop = FALSE)
 #head(y)
#check if wy is "variable" or not: argument wy_rg
##this lines are creating a wy with no interactions within regimes
  if(wy_rg){
    wy <- Matrix(0, nrow = n, ncol = sv)
for(i in 1:sv) wy[,i] <- (Ws %*% (as.matrix(y)*rgm[,i]))*rgm[,i]
colnames(wy) <- paste(paste("W_",colnames(y)),rep(1:sv), sep = "_")
  }
else{
  wy <- Ws %*% as.matrix(y)
  colnames(wy) <- paste("W_", colnames(y), sep ="")
}

#head(wy)

 ###############################
  ###############################
  ##### First two parts of formula: x fixed and x variable
  ###############################
  ###############################

   ## extract X fixed
  Xf <- model.matrix(F1, data = mf, rhs = 1, drop = FALSE)
  namesxf <- colnames(Xf)
  ## extract X variable
  Xv <- model.matrix(F1, data = mf, rhs = 2, drop = FALSE)
  namesxv <- colnames(Xv)


  if(any(namesxf == "(Intercept)") && any(namesxv == "(Intercept)"))
    stop("(Intercept) cannot  be specified as fixed and variable regressor at the same time!")
  if(!("(Intercept)" %in% c(namesxf, namesxv)))
    warning("The model has been specified without an intercept")

  k2 <- dim(Xv)[2]
  totc <- sv*k2
  XV <- matrix(0, ncol = totc, nrow = n)
  seq_1 <- seq(1, totc, k2)
  seq_2 <- seq(k2, totc,  k2)
  for(i in 1:sv) XV[,seq_1[i]:seq_2[i]] <- Xv * rgm[,i]
  namesxV <- paste(namesxv,rep(1:sv,each = k2), sep = "_")
  colnames(XV) <- namesxV
 # head(Xf)
  #head(XV)


###############################
  ###############################
  ##### Formula 3  wx (and instruments)
  ###############################
  ###############################

  wx <- model.matrix(F1, data = mf, rhs = 3, drop = FALSE)


  if(any(colnames(wx) == "(Intercept)")) wx <- wx[,-which(colnames(wx) == "(Intercept)")]
  nameswx <- colnames(wx)
 # head(wx)

  xfd  <- wx[,(nameswx %in% namesxf), drop = FALSE]
  xvd  <- wx[,(nameswx %in%  namesxv), drop = FALSE]
  namesxfd <- colnames(xfd)
  namesxvd <- colnames(xvd)

  #head(xfd)
  #head(xvd)
  nameswxf <- paste("W_",namesxfd, sep="")
#   nameswwxf <- paste("WW_",namesxfd, sep="")
# ##just for the instruments
#   nameswwwxf <- paste("WWW_",namesxfd, sep="")

  ### if x is fixed wx is fixed
  Wxf <- Ws %*% xfd
  colnames(Wxf) <- nameswxf
  # WWxf <- Ws %*% Wxf
  # colnames(WWxf) <- nameswwxf
  # WWWxf <- Ws %*% WWxf
  # colnames(WWWxf) <- nameswwwxf
#head(Xf)
#head(Wxf)
#head(WWxf)
#head(WWWxf)

### if x varies wx varies
  #nameswxv <- paste(paste("W_",colnames(wxv), sep=""), rep(1:sv,each = k3), sep = "_")

namesxvD <- paste(namesxvd,rep(1:sv,each = ncol(xvd)), sep = "_")
xvD  <- XV[, which(namesxV %in% namesxvD), drop = FALSE]

#head(xvD)
#head(XV)

####take into account that xvD is a matrix
  WxvD <- matrix(0, ncol = ncol(xvD), nrow = n)
  seq_1 <- seq(1, ncol(xvD), sv)
  seq_2 <- seq(sv, ncol(xvD),  sv)

  for (i in 1: (ncol(xvD)/sv)) WxvD[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (xvD[,seq_1[i]:seq_2[i]]*rgm[,i]))*rgm[,i])
  nameswxv <- paste("W_",colnames(xvD), sep="")
  colnames(WxvD) <- nameswxv

 # head(xvD)
#  head(WxvD)



  Zmat <- cbind(Xf, XV, Wxf, WxvD, wy)
#head(Zmat)



###############################
###############################
##### parts 4 and 5 of formula: instrument fixed and variable
###############################
###############################



  Zf <- model.matrix(F1, data = mf, rhs = 4, drop = FALSE)
  nameszf <- colnames(Zf)
  Zv <- model.matrix(F1, data = mf, rhs = 5, drop = FALSE)
  nameszv <- colnames(Zv)
  namesH <- c(nameszf, nameszv)


  #head(Zf)
  #head(Zv)
  if(any(nameszf == "(Intercept)") && any(nameszv == "(Intercept)"))
    stop("(Intercept) cannot  be specified as fixed and variable instruments at the same time!")

  if(any(namesxf == "(Intercept)") && any(nameszv == "(Intercept)"))
    stop("(Intercept) in instruments should not vary")

  if(any(namesxv == "(Intercept)") && any(nameszf == "(Intercept)"))
    stop("(Intercept) in instruments should vary")

  if(!("(Intercept)" %in% c(nameszf, nameszv)))
    warning("The model has been specified without an intercept but the instruments include the intercept")

  k3 <- dim(Zv)[2]
  totz <- sv*k3
  ZV <- matrix(0, ncol = totz, nrow = n)
  seq_1 <- seq(1, totz, k3)
  seq_2 <- seq(k3, totz,  k3)
  for(i in 1:sv) ZV[,seq_1[i]:seq_2[i]] <- Zv * rgm[,i]
  nameszv <- paste(nameszv,rep(1:sv,each = k3), sep = "_")
  colnames(ZV) <- nameszv

  #head(Zf)
  #head(ZV)

  #### check for the endogenous for identification and printing
  end.f <- Xf[, !(colnames(Xf) %in% colnames(Zf)), drop = FALSE]
  end.v <- Xv[, !(colnames(Xv) %in% colnames(Zv)), drop = FALSE]
  colnames.end.f <- colnames(end.f)
  colnames.end.v <- colnames(end.v)
  namesinstr.F <- colnames(Zf[ , !(colnames(Zf) %in% colnames(Xf)), drop = FALSE])
  namesinstr.V <- colnames(Zv[ , !(colnames(Zv) %in% colnames(Xv)), drop = FALSE])

  namesinstr.VV <-  paste(namesinstr.V, rep(1:sv,each = length(namesinstr.V)) , sep = "_")
  nameInsts <- c(namesinstr.F, namesinstr.V)
  nameInst <- c(namesinstr.F, namesinstr.VV)
  #print(nameInst)
 # head(end.f)
  #head(end.v)
  end.fl <- xfd[, (colnames(xfd) %in% colnames(end.f)), drop = FALSE]
  end.vl <- xvd[, (colnames(xvd) %in% colnames(end.v)), drop = FALSE]
  if(!is.null(colnames(end.fl))) colnames.end.fl <- paste("W_",colnames(end.fl), sep ="")
  if(!is.null(colnames(end.vl))) colnames.end.vl <-paste("W_",paste(colnames(end.vl), rep(1:sv,each = length(colnames(end.vl))), sep = "_"), sep = "")
#head(end.fl)
#head(end.vl)
  ### x fixed and variables for the instruments
  x.f <- Xf[, (colnames(Xf) %in% colnames(Zf)), drop = FALSE]
  x.v <- Xv[, (colnames(Xv) %in% colnames(Zv)), drop = FALSE]

  if(any(colnames(x.f) == "(Intercept)")) x.f <- x.f[,-which(colnames(x.f) == "(Intercept)")]
  if(any(colnames(x.v) == "(Intercept)")) x.v <- x.v[,-which(colnames(x.v) == "(Intercept)")]
 # head(x.f)
#  head(x.v)
  ##need to build the instruments
  namesx.f <- colnames(x.f)
  namesx.v <- colnames(x.v)
  nameswx.f <- paste("W_",namesx.f, sep="")
  nameswwx.f <- paste("WW_",namesx.f, sep="")

  Wx.f <- Ws %*% x.f
  colnames(Wx.f) <- nameswx.f
  WWx.f <- Ws %*% Wx.f
  colnames(WWx.f) <- nameswwx.f
  #only for the durbin variables third power
   WWWx.f <- Ws %*% WWx.f[,which(nameswx %in% namesxfd), drop = F]
   colnames(WWWx.f) <- paste("WWW_",namesx.f, sep="")

   x.f.inst <- Xf[,-which(namesxf %in% colnames.end.f)]

   Hx.fne <- cbind(x.f.inst, Wx.f, WWx.f, WWWx.f)
 #  head(Hx.fne)
   ### now the x.v
   namesx.V <-  paste(namesx.v, rep(1:sv, each = length(namesx.v)), sep = "_")

   x.V <- XV[,which(namesxV %in% namesx.V)]
  # head(x.V)

   Wx.V <- matrix(0, ncol = ncol(x.V), nrow = n)
   seq_1 <- seq(1, ncol(x.V), sv)
   seq_2 <- seq(sv, ncol(x.V),  sv)

   for (i in 1: (ncol(x.V)/sv)) Wx.V[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (x.V[,seq_1[i]:seq_2[i]]*rgm[,i]))*rgm[,i])
   nameswx.V <- paste("W_",colnames(x.V), sep="")
   colnames(Wx.V) <- nameswx.V
    #head(Wx.V)
    WWx.V <- matrix(0, ncol = ncol(x.V), nrow = n)
    for (i in 1: (ncol(x.V)/sv)) WWx.V[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (Wx.V[,seq_1[i]:seq_2[i]]*rgm[,i]))*rgm[,i])
    nameswwx.V <- paste("WW_",colnames(x.V), sep="")
    colnames(WWx.V) <- nameswwx.V

    WWx.Vfd <- WWx.V[,which(namesx.V %in% colnames(xvD))]
    seq_1 <- seq(1, ncol(WWx.Vfd), sv)
    seq_2 <- seq(sv, ncol(WWx.Vfd),  sv)
    WWWx.V <- matrix(0, ncol = ncol(WWx.Vfd), nrow = n)
    for (i in 1: (ncol(WWx.Vfd)/sv)) WWWx.V[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (WWx.Vfd[,seq_1[i]:seq_2[i]]*rgm[,i]))*rgm[,i])
    nameswwwx.V <- paste("WWW_",colnames(WWx.Vfd), sep="")
    colnames(WWWx.V) <- nameswwwx.V

    colnames.end.V <- paste(colnames.end.v, rep(1:sv, each = length(colnames.end.v)), sep = "_")
    x.v.inst <- XV[,-which(namesxV %in% colnames.end.V)]
    #head(x.v.inst)

    Hx.vne <- cbind(x.v.inst, Wx.V, WWx.V, WWWx.V)
  ###############################
  ###############################
  ##### Formula 6  winst
  ###############################
  ###############################


  wh <- model.matrix(F1, data = mf, rhs = 6, drop = FALSE)
 colnwh <- colnames(wh)
 if(!any(colnames(wh) %in% nameInsts))
   stop("This part of the formula is only for external instruments")

 if(any(colnwh == "(Intercept)")){
    wh <- as.matrix(wh[,-which(colnwh  == "(Intercept)")])
    colnames(wh) <- colnwh[-which(colnwh == "(Intercept)")]
  }

 whf  <- wh[,(colnames(wh) %in% colnames(Zf)), drop = FALSE]
 whv  <- wh[,(colnames(wh) %in%  colnames(Zv)), drop = FALSE]
 nameswhf <- colnames(whf)
 nameswhv <- colnames(whv)


 winsf <- as.matrix(Ws %*% whf)

 if(dim(winsf)[2]!=0){
 nameswhf <- paste("W_",nameswhf, sep="")
 colnames(winsf) <- nameswhf
}

 nameswhV <- paste(nameswhv,rep(1:sv,each = ncol(whv)), sep = "_")
 hvD  <- ZV[, which(nameszv %in% nameswhV), drop = FALSE]



 WhvD <- matrix(0, ncol = ncol(hvD), nrow = n)
 seq_1 <- seq(1, ncol(hvD), sv)
 seq_2 <- seq(sv, ncol(hvD),  sv)

 for (i in 1: (ncol(hvD)/sv)) WhvD[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (hvD[,seq_1[i]:seq_2[i]]*rgm[,i]))*rgm[,i])
 nameswhv <- paste("W_",colnames(hvD), sep="")
 colnames(WhvD) <- nameswhv

 #head(WhvD)

 Hmat <- cbind(Hx.fne, Hx.vne, winsf, WhvD)

 Hmat <- Hmat[, qr(Hmat)$pivot[seq_len(qr(Hmat)$rank)]]


## if instruments check identification


  colinst <- c(nameInst, nameswhf, nameswhv)

  colnames.end.V <-  paste(colnames.end.v, rep(1:sv,each = length(colnames(end.vl))) , sep = "_")
  col.end.f.l <- colnames(end.fl)
  colnames.end.V.l <-  paste("W_",paste(colnames(end.vl), rep(1:sv,each = length(colnames(end.vl))) , sep = "_"), sep="")

  if(length(colinst) < length(c(colnames.end.f, colnames.end.V, col.end.f.l, colnames.end.V.l)))
    stop("Not enough instruments specified: the model is not identified")

  ret <- list(y = y, Hmat = Hmat, Zmat = Zmat, endog = list(colnames.end.f, colnames.end.V, col.end.f.l, colnames.end.V.l),
              instrum = colinst)
  return(ret)
}


spatial.ivreg.regimes <-function(y, Zmat, Hmat, het){
  df <- nrow(Zmat) - ncol(Zmat)
  HH <- crossprod(Hmat,Hmat)
  Hye <- crossprod(Hmat, Zmat)
  bz <- solve(HH,Hye)
  Zp <- Hmat %*% bz
  ZpZp <- crossprod(Zp)
  ZpZpi <- solve(ZpZp)
  Zpy <- crossprod(Zp,y)
  delta <- crossprod(ZpZpi,Zpy)
  yp <- Zmat %*% delta
  e <- y - yp
  #	print(dim(e))

    if(het)	{

      s2 <- crossprod(e) /df
      omega <- as.numeric(e^2)
      ZoZ <- crossprod(Zp, (Zp * omega))
      vardelta <- ZpZpi %*% ZoZ %*% ZpZpi

    }
    else{
      s2 <- crossprod(e) / df
      vardelta <- ZpZpi * as.numeric(s2)
    }

  result <- list(coefficients = delta, var = vardelta)
  return(result)
}

