
lag_regimes <- function(formula, data, listw, rgv,
                   het, cl, wy_rg){




  intro <- iv.lag.data.prep.regimes(formula, data = data, listw = listw,
                                    wy_rg = wy_rg, rgv = rgv, weps_rg = FALSE)


  y        <- intro[[1]]
  Hmat     <- intro[[2]]
  Zmat     <- intro[[3]]
  colnames.end<- intro[[4]]
  colnames.instr <- intro[[5]]
  colinstr <- intro[[8]]
  nameswx  <- intro[[9]]
  res <- spatial.ivreg.regimes(as.matrix(y), as.matrix(Zmat), as.matrix(Hmat), het)
  res <- list(res, cl, colnames.end,  colnames.instr, colinstr, nameswx)
  class(res) <- c("spregimes","lag_regimes")
  return(res)

}




iv.lag.data.prep.regimes <- function(formula, data, rgv, listw, wy_rg, weps_rg){


  #define the split variable
  if(is.null(rgv)) stop("regimes variable not specified")
  if(!inherits(rgv, "formula")) stop("regimes variable has to be a formula")

  n                <- dim(data)[1]
  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- as.numeric(unique(splitvar))
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))

  l.split <- list(n, splitvar, sv, svm, rgm)
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

  if(weps_rg && attributes(F1)$rhs[[1]]!=0) stop("If the spatial process vary by
                                                                      regimes, all of the regressors should vary by regimes")

  mf <- model.frame(F1, data = data)
  y <- model.part(F1, data = mf, lhs = 1, drop = FALSE)

#check if wy is "variable" or not: argument wy_rg

  if(wy_rg){
    wy <- Matrix(0, nrow = n, ncol = sv)
    for(i in 1:sv) wy[,i] <- Ws %*% (as.matrix(y)*rgm[,i])
    colnames(wy) <- paste("W_",paste(colnames(y),rep(1:sv), sep = "_"), sep="")
  }
else{
    wy <- Ws %*% as.matrix(y)
    colnames(wy) <- paste("W_", colnames(y), sep ="")
}


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
  namesxv <-  colnames(Xv)

  if(any(namesxf == "(Intercept)") && any(namesxv == "(Intercept)"))
    stop("(Intercept) cannot  be specified as fixed and variable regressor at the same time!")
  if(!("(Intercept)" %in% c(namesxf, namesxv)))
    warning("The model has been specified without an intercept")

if(dim(Xv)[2] != 0){
  k2 <- dim(Xv)[2]
  totc <- sv*k2
  XV <- matrix(0, ncol = totc, nrow = n)
  seq_1 <- seq(1, totc, k2)
  seq_2 <- seq(k2, totc,  k2)
  for(i in 1:sv) XV[,seq_1[i]:seq_2[i]] <- Xv * rgm[,i]
  namesxV <- paste(namesxv,rep(1:sv,each = k2), sep = "_")
  colnames(XV) <- namesxV
}
else{
  XV <- matrix(0, ncol = 0, nrow = n)
}
  ###############################
  ###############################
  ##### Formula 3  wx
  ###############################
  ###############################
  # variables for Durbin

  wx <- model.matrix(F1, data = mf, rhs = 3, drop = FALSE)
  nameswx <-  colnames(wx)

  if(any(nameswx == "(Intercept)")) wx <- wx[,-which(colnames(wx) == "(Intercept)")]
  if(any(nameswx == "(Intercept)")) nameswx <- nameswx[-(which(nameswx == "(Intercept)"))]
  wx <- as.matrix(wx)
  colnames(wx) <- nameswx

#check if Durbin are fixed or variable
  xfd  <- wx[,(nameswx %in% namesxf), drop = FALSE]
  xvd  <- wx[,(nameswx %in%  namesxv), drop = FALSE]
  namesxfd <- colnames(xfd)
  namesxvd <- colnames(xvd)
  namesd <- c(namesxfd, namesxvd)
  ### if x is fixed wx is fixed
 if(dim(xfd)[2] !=0){
  Wxf <- Ws %*% xfd
  colnames(Wxf) <- paste("W_",namesxfd, sep="")
  }
  else {
    Wxf <- matrix(nrow = n, ncol = 0)
    colnames(Wxf) <- NULL
  }

  ### if x varies wx varies
if(dim(xvd)[2] != 0){
  namesxvD <- paste(namesxvd,rep(1:sv,each = ncol(xvd)), sep = "_")
  xvD  <- XV[, which(namesxV %in% namesxvD), drop = FALSE]
  WxvD <- matrix(0, ncol = ncol(xvD), nrow = n)
  seq_1 <- seq(1, ncol(xvD), ncol(xvD)/sv)
  seq_2 <- seq(ncol(xvD)/sv, ncol(xvD),  ncol(xvD)/sv)
for (i in 1: sv) WxvD[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (xvD[,seq_1[i]:seq_2[i]])))
  nameswxv <- paste("W_",colnames(xvD), sep="")
  colnames(WxvD) <- nameswxv
}
  else{
    namesxvD <- NULL
    xvD <- matrix(nrow = n, ncol = 0)
    WxvD <- matrix(nrow = n, ncol = 0)
  }

  Zmat <- cbind(Xf, XV, Wxf, WxvD, wy)


###############################
###############################
##### parts 4 and 5 of formula: instrument, endogenous and x no.end
###############################
###############################


#### Instruments:
  ### 1) endog in X but not in Z
  ### 2) external instr not in x but only in z
  ### 3) esog in both

  Zf <- model.matrix(F1, data = mf, rhs = 4, drop = FALSE)
  nameszf <-  colnames(Zf)
  Zv <- model.matrix(F1, data = mf, rhs = 5, drop = FALSE)
  nameszv <-  colnames(Zv)
  namesH <- c(nameszf, nameszv)

  if(any(nameszf == "(Intercept)") && any(nameszv == "(Intercept)"))
    stop("(Intercept) cannot  be specified as fixed and variable instruments at the same time!")

  if(any(namesxf == "(Intercept)") && any(nameszv == "(Intercept)"))
    stop("(Intercept) in instruments should not vary")

  if(any(namesxv == "(Intercept)") && any(nameszf == "(Intercept)"))
    stop("(Intercept) in instruments should vary")

  if(!("(Intercept)" %in% c(nameszf, nameszv)))
    warning("The model has been specified without an intercept but the instruments include the intercept")

  if(dim(Zv)[2] != 0){

    k3 <- dim(Zv)[2]
    totz <- sv*k3
    ZV <- matrix(0, ncol = totz, nrow = n)
    seq_1 <- seq(1, totz, k3)
    seq_2 <- seq(k3, totz,  k3)
    for(i in 1:sv) ZV[,seq_1[i]:seq_2[i]] <- Zv * rgm[,i]
    nameszv <- paste(nameszv,rep(1:sv,each = k3), sep = "_")
    colnames(ZV) <- nameszv
  }
  else{
    ZV <- matrix(nrow = n, ncol = 0)
    nameszv <- NULL
  }

  ### 1) endog in X but not in Z (only needed for summary)
  end.f <- Xf[, !(colnames(Xf) %in% colnames(Zf)), drop = FALSE]
  end.v <- Xv[, !(colnames(Xv) %in% colnames(Zv)), drop = FALSE]
  colnames.end.f <- colnames(end.f)
  colnames.end.v <- colnames(end.v)
  #endogenous lagged (only needed for summary)
  end.fl <- xfd[, (colnames(xfd) %in% colnames(end.f)), drop = FALSE]
  end.vl <- xvd[, (colnames(xvd) %in% colnames(end.v)), drop = FALSE]
  col.end.f.l <- colnames(end.fl)
  if (dim(end.v)[2] != 0) colnames.end.V <- paste(colnames.end.v, rep(1:sv, each = length(colnames(end.vl))) , sep = "_")
  else colnames.end.V <- NULL
  if(!is.null(colnames(end.fl))) colnames.end.fl <- paste("W_",colnames(end.fl), sep ="")
  else colnames.end.fl <- NULL
  if(!is.null(colnames(end.vl))) colnames.end.vl <- paste("W_", paste(colnames(end.vl), rep(1:sv,each = length(colnames(end.vl))), sep = "_"), sep = "")
  else colnames.end.vl <- NULL


  ### 2) external instr not in x but only in z
  instr.F <- Zf[ , !(colnames(Zf) %in% colnames(Xf)), drop = FALSE]
  instr.V <- Zv[ , !(colnames(Zv) %in% colnames(Xv)), drop = FALSE]
  if(dim(instr.V)[2] != 0){
    k3 <- dim(instr.V)[2]
    totz <- sv*k3
    instr.VV <- matrix(0, ncol = totz, nrow = n)
    seq_1 <- seq(1, totz, k3)
    seq_2 <- seq(k3, totz,  k3)
    for(i in 1:sv) instr.VV[,seq_1[i]:seq_2[i]] <- instr.V * rgm[,i]
    colnames(instr.VV) <- paste(colnames(instr.V), rep(1:sv,each = k3), sep = "_")
  }
  else{
    instr.VV <- matrix(nrow = n, ncol = 0)
    colnames(instr.VV) <- NULL
  }
  namesinstr.F <- colnames(Zf[ , !(colnames(Zf) %in% colnames(Xf)), drop = FALSE])
  namesinstr.V <- colnames(Zv[ , !(colnames(Zv) %in% colnames(Xv)), drop = FALSE])
  if(length(namesinstr.V) !=0) namesinstr.VV <-  paste(namesinstr.V, rep(1:sv,each = length(namesinstr.V)) , sep = "_")
  else namesinstr.VV <- NULL
  ##needed  for printing and check identification
  nameInsts <- c(namesinstr.F, namesinstr.V)
  nameInst <- c(namesinstr.F, namesinstr.VV)


  ### 3) esog in both
  x.ff <- Xf[, (colnames(Xf) %in% colnames(Zf)), drop = FALSE]
  x.vv <- Xv[, (colnames(Xv) %in% colnames(Zv)), drop = FALSE]
 if(any(colnames(x.ff) == "(Intercept)")) x.f <- x.ff[,-which(colnames(x.ff) == "(Intercept)")]
  else x.f <- x.ff
  if(any(colnames(x.vv) == "(Intercept)")) x.v <- x.vv[,-which(colnames(x.vv) == "(Intercept)")]
  else x.v <- x.vv
  if(dim(x.f)[2] != 0){
    namesx.f <- colnames(x.f)
    nameswx.f <- paste("W_",namesx.f, sep="")
    nameswwx.f <- paste("WW_",namesx.f, sep="")
    Wx.f <- Ws %*% x.f
    colnames(Wx.f) <- nameswx.f
    WWx.f <- Ws %*% Wx.f
    colnames(WWx.f) <- nameswwx.f
  if(dim(xfd)[2] !=0){
   WWWx.f <- Ws %*% WWx.f[,which(nameswx %in% namesxfd), drop = F]
   WWWx.f <- as.matrix(WWWx.f)
   colnames(WWWx.f) <- paste("WWW_", namesxfd, sep="")
  }
  else WWWx.f <- matrix(nrow = n, ncol = 0)
  }
  else {
    Wx.f <-  WWx.f <- WWWx.f <- matrix(nrow = n, ncol = 0)
    namesx.f <-  namesx.v <- nameswx.f <- nameswwx.f <- NULL
  }
Hx.fne <- cbind(x.ff, Wx.f, WWx.f, WWWx.f, instr.F)

 if(dim(x.v)[2] != 0 ){
     namesx.vv <- colnames(x.vv)
     namesx.VV <-  paste(namesx.vv, rep(1:sv, each = length(namesx.vv)), sep = "_")
     x.VV <- XV[,which(namesxV %in% namesx.VV), drop = F]
     namesx.v <- colnames(x.v)
     namesx.V <-  paste(namesx.v, rep(1:sv, each = length(namesx.v)), sep = "_")
     x.V <- XV[,which(namesxV %in% namesx.V), drop = F]
     Wx.V <- matrix(0, ncol = ncol(x.V), nrow = n)
     seq_1 <- seq(1, ncol(x.V), (ncol(x.V)/sv))
     seq_2 <- seq(ncol(x.V)/sv, ncol(x.V),ncol(x.V)/sv)
     for(i in 1:sv)  Wx.V[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (x.V[,seq_1[i]:seq_2[i]])))
     nameswx.V <- paste("W_",colnames(x.V), sep="")
     colnames(Wx.V) <- nameswx.V
     WWx.V <- matrix(0, ncol = ncol(x.V), nrow = n)
     for (i in 1: sv) WWx.V[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (Wx.V[,seq_1[i]:seq_2[i]])))
     nameswwx.V <- paste("WW_",colnames(x.V), sep="")
     colnames(WWx.V) <- nameswwx.V
     WWx.Vfd <- WWx.V[,which(namesx.V %in% colnames(xvD))]
       if(dim(WWx.Vfd)[2] != 0){
            seq_1 <- seq(1, ncol(WWx.Vfd),ncol(WWx.Vfd)/ sv)
            seq_2 <- seq(ncol(WWx.Vfd)/ sv, ncol(WWx.Vfd),  ncol(WWx.Vfd)/ sv)
            WWWx.V <- matrix(0, ncol = ncol(WWx.Vfd), nrow = n)
              for (i in 1:sv) WWWx.V[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (WWx.Vfd[,seq_1[i]:seq_2[i]])))
            nameswwwx.V <- paste("W",colnames(WWx.Vfd), sep="")
            colnames(WWWx.V) <- nameswwwx.V
     }
    else WWWx.V <- matrix(nrow = n, ncol = 0)
    if(length(colnames.end.V) != 0) colnames.end.V <- paste(colnames.end.v, rep(1:sv, each = length(colnames.end.v)), sep = "_")
    else colnames.end.V <- NULL
   }
   else{
     x.VV <- Wx.V <- WWx.V <-  WWWx.V <- matrix(nrow = n, ncol = 0)
     colnames.end.V <- paste(colnames.end.v, rep(1:sv, each = length(colnames.end.v)), sep = "_")
   }
Hx.vne <- cbind(x.VV, Wx.V, WWx.V, WWWx.V, instr.VV)

  ###############################
  ###############################
  ##### Formula 6  winst
  ###############################
  ###############################


  wh <- model.matrix(F1, data = mf, rhs = 6, drop = FALSE)

     if(dim(wh)[2] != 0){
          colnwh <- colnames(wh)
               if(!any(colnwh %in% nameInsts))
                      stop("This part of the formula is only for external instruments")

            if(any(colnwh == "(Intercept)")){
                  wh <- as.matrix(wh[,-which(colnwh  == "(Intercept)")])
                      colnames(wh) <- colnwh[-which(colnwh == "(Intercept)")]
                      }

 whf  <- wh[,(colnames(wh) %in% colnames(Zf)), drop = FALSE]
 whv  <- wh[,(colnames(wh) %in%  colnames(Zv)), drop = FALSE]


 if(dim(whf)[2] != 0){
   nameswhf <- colnames(whf)
   winsf <- as.matrix(Ws %*% whf)
   nameswhf <- paste("W_",nameswhf, sep="")
   colnames(winsf) <- nameswhf
}
else{
  winsf <- matrix(nrow = n, ncol = 0 )
  nameswhf <- NULL
}

if(dim(whv)[2] != 0){
 nameswhv <- colnames(whv)
 nameswhV <- paste(nameswhv,rep(1:sv,each = ncol(whv)), sep = "_")
 hvD  <- ZV[, which(nameszv %in% nameswhV), drop = FALSE]
 WhvD <- matrix(0, ncol = ncol(hvD), nrow = n)
 seq_1 <- seq(1, ncol(hvD), sv)
 seq_2 <- seq(sv, ncol(hvD),  sv)
 for (i in 1: sv) WhvD[,seq_1[i]:seq_2[i]] <-  as.matrix((Ws %*% (hvD[,seq_1[i]:seq_2[i]])))
 nameswhv <- paste("W_",colnames(hvD), sep="")
 colnames(WhvD) <- nameswhv
    }
    else {
      WhvD <- matrix(nrow = n, ncol = 0)
      nameswhv <- NULL
    }
}
 else{
   winsf <- matrix(nrow = n, ncol = 0 )
   nameswhf <- NULL
   WhvD <- matrix(nrow = n, ncol = 0)
   nameswhv <- NULL
 }

 Hmat <- cbind(Hx.fne, Hx.vne, winsf, WhvD)
 Hmat <- Hmat[, qr(Hmat)$pivot[seq_len(qr(Hmat)$rank)]]

colinstr <- c(nameInst, nameswhf, nameswhv)


if(!is.null(colinstr)) {

    if(length(colinstr) < length(c(colnames.end.f, colnames.end.V, col.end.f.l, colnames.end.vl)))
    stop("Not enough instruments specified: the model is not identified")

  if(!is.null(namesd)) colinst <- c("X", "WX","WWX", "WWWX", colinstr)
  else colinst <- c("X", "WX","WWX",  colinstr)
}
  else {
    if(!is.null(namesd)) colinst <- c("X", "WX","WWX","WWWX")
    else colinst <- c("X", "WX","WWX")
    }

 ret <- list(y = y, Hmat = Hmat, Zmat = Zmat, endog = list(colnames(wy), colnames.end.f, colnames.end.V, col.end.f.l, colnames.end.vl),
              instrum = colinst, l.split = l.split, ws = Ws, colinstr = colinstr, nameswx = nameswx)
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
  yp <- array(yp, dim = c(length(yp),1),
              dimnames = list(seq(1, length(yp)), ""))
  e <- y - yp
  e <- array(e, dim = c(length(e),1),
             dimnames = list(seq(1, length(e)), ""))
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

  result <- list(coefficients = delta, var = vardelta,
                 residuals = e,
                 X = Zmat, y = y, yp = yp)
  return(result)
}


