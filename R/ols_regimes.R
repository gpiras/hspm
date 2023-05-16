


ols_regimes <- function(formula, data, listw, rgv,
                        het, cl){




  intro <- ols.data.prep.regimes(formula, data = data, listw = listw, rgv = rgv,
                                 weps_rg = FALSE, model = "ols")


  y        <- intro[[1]]
  Hmat     <- intro[[2]]
  Zmat     <- intro[[3]]
  colnames.end   <- intro[[4]]
  colnames.instr <- intro[[5]]
  dur       <- intro[[9]]
  ct <- intro[[6]][[6]]

  res <- spatial.ivreg.regimes(as.matrix(y), as.matrix(Zmat), as.matrix(Hmat), het)

#print(colnames.end)
#print(colnames.instr)
  res <- Matchgroups(res, ct)
  colnames.end <- Matchnames(colnames.end, ct)
  colnames.instr <- Matchnames(colnames.instr, ct)
  res <- list(res, cl, colnames.end,  colnames.instr, dur)
  class(res) <- c("spregimes", "ols_regimes")


  return(res)

}


ols.data.prep.regimes <- function(formula, data, rgv, listw,
                                  weps_rg, model){


  #define the split variable
  if(is.null(rgv)) stop("regimes variable not specified")
  if(!inherits(rgv, "formula")) stop("regimes variable has to be a formula")

  n                <- dim(data)[1]
  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- as.numeric(unique(splitvar))
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))
  mt                <- cbind(1:sv, as.numeric(unique(splitvar)))
  ct                <- as.numeric(mt[order(mt[,2], decreasing = F),1])
  ct                <- cbind(ct, svm)

  l.split <- list(n, splitvar, sv, svm, rgm, ct)

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
    namesxV <- paste(namesxv, rep(1:sv, each = k2), sep = "_")
    colnames(XV) <- namesxV
  }

  ###############################
  ###############################
  ##### Formula 3  wx
  ###############################
  ###############################
  # variables for Durbin

  wx <- model.matrix(F1, data = mf, rhs = 3, drop = FALSE)
  nmwx <- colnames(wx)
  if(any(nmwx == "(Intercept)")) {
    wx <- as.matrix(wx[,-which(colnames(wx) == "(Intercept)")])
    colnames(wx) <- nmwx[-which(nmwx == "(Intercept)")]
  }
  if(ncol(wx) != 0) dur <- TRUE
  else dur <- FALSE
   nameswx <-  colnames(wx)
  #check if Durbin are fixed or variable
  xfd  <- wx[,(nameswx %in% namesxf), drop = FALSE]
  xvd  <- wx[,(nameswx %in%  namesxv), drop = FALSE]
  namesxfd <- colnames(xfd)
  namesxvd <- colnames(xvd)

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
    namesxvD <- paste(namesxvd,rep(1:sv, each = ncol(xvd)), sep = "_")
    xvD  <- XV[, which(namesxV %in% namesxvD), drop = FALSE]
    WxvD <- matrix(0, ncol = ncol(xvD), nrow = n)
    seq_1 <- seq(1, ncol(xvD), ncol(xvD)/sv)
    seq_2 <- seq(ncol(xvD)/sv, ncol(xvD),  ncol(xvD)/sv)
    for (i in 1: sv) WxvD[,seq_1[i]:seq_2[i]] <-  as.matrix(Ws) %*% xvD[,seq_1[i]:seq_2[i]]
    nameswxv <- paste("W_",colnames(xvD), sep="")
    colnames(WxvD) <- nameswxv
  }
  else{
    namesxvD <- NULL
    xvD <- matrix(nrow = n, ncol = 0)
    WxvD <- matrix(nrow = n, ncol = 0)
  }

Zmat <- cbind(Xf, XV, Wxf, WxvD)

  ###############################
  ###############################
  ##### parts 4 and 5 of formula: instrument, endogenous and x no.end
  ###############################
  ###############################


  #### Instruments:
  ### 1) endog in X but not in Z
  ### 2) external instr not in x but only in z
  ### 3) esog in both: this model doesn't need this part

  Zf <- model.matrix(F1, data = mf, rhs = 4, drop = FALSE)
  nameszf <-  colnames(Zf)
  Zv <- model.matrix(F1, data = mf, rhs = 5, drop = FALSE)
  nameszv <-  colnames(Zv)
  namesH <- c(nameszf, nameszv)

if(!is.null(namesH)){
  if(any(nameszf == "(Intercept)") && any(nameszv == "(Intercept)"))
    stop("(Intercept) cannot  be specified as fixed and variable instruments at the same time!")

  if(any(namesxf == "(Intercept)") && any(nameszv == "(Intercept)"))
    stop("(Intercept) in instruments should not vary")

  if(any(namesxv == "(Intercept)") && any(nameszf == "(Intercept)"))
    stop("(Intercept) in instruments should vary")

  if(!("(Intercept)" %in% c(nameszf, nameszv)))
    warning("The model has been specified without an intercept
            but the instruments include the intercept")

  if(dim(Zv)[2] != 0){

    k3 <- dim(Zv)[2]
    totz <- sv*k3
    ZV <- matrix(0, ncol = totz, nrow = n)
    seq_1 <- seq(1, totz, k3)
    seq_2 <- seq(k3, totz,  k3)
    for(i in 1:sv) ZV[,seq_1[i]:seq_2[i]] <- Zv * rgm[,i]
    nameszv <- paste(nameszv,rep(1:sv, each = k3), sep = "_")
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
  if (dim(end.v)[2] != 0) colnames.end.V <- paste(colnames.end.v,
                                                  rep(1:sv, each = length(colnames(end.vl))) , sep = "_")
  else colnames.end.V <- NULL
  if(!is.null(colnames(end.fl))) colnames.end.fl <- paste("W_",colnames(end.fl), sep ="")
  else colnames.end.fl <- NULL
  if(!is.null(colnames(end.vl))) colnames.end.vl <- paste("W_", paste(colnames(end.vl),
                                                                      rep(1:sv, each = length(colnames(end.vl))), sep = "_"), sep = "")
  else colnames.end.vl <- NULL


  ### 2) external instr not in x but only in z
  # add the external instr
  instr.F <- Zf[ , !(colnames(Zf) %in% colnames(Xf)), drop = FALSE]
  instr.V <- Zv[ , !(colnames(Zv) %in% colnames(Xv)), drop = FALSE]

  if(dim(instr.V)[2] != 0){
    k3 <- dim(instr.V)[2]
    totz <- sv*k3
    instr.VV <- matrix(0, ncol = totz, nrow = n)
    seq_1 <- seq(1, totz, k3)
    seq_2 <- seq(k3, totz,  k3)
    for(i in 1:sv) instr.VV[,seq_1[i]:seq_2[i]] <- instr.V * rgm[,i]
    colnames(instr.VV) <- paste(colnames(instr.V), rep(1:sv, each = k3), sep = "_")
  }
  else{
    instr.VV <- matrix(nrow = n, ncol = 0)
    colnames(instr.VV) <- NULL
  }
  namesinstr.F <- colnames(Zf[ , !(colnames(Zf) %in% colnames(Xf)), drop = FALSE])
  namesinstr.V <- colnames(Zv[ , !(colnames(Zv) %in% colnames(Xv)), drop = FALSE])
  if(length(namesinstr.V) !=0) namesinstr.VV <-  paste(namesinstr.V,
                                                       rep(1:sv, each = length(namesinstr.V)) , sep = "_")
  else namesinstr.VV <- NULL
  ##need to for pinting and check identification
  nameInsts <- c(namesinstr.F, namesinstr.V)
  nameInst <- c(namesinstr.F, namesinstr.VV)

  ### 3) esog in both
  x.f <- Xf[, (colnames(Xf) %in% colnames(Zf)), drop = FALSE]
  x.v <- Xv[, (colnames(Xv) %in% colnames(Zv)), drop = FALSE]
Hx.fne <- cbind(x.f, Wxf,instr.F)
  if(dim(x.v)[2] != 0 ){
    namesx.v <- colnames(x.v)
    namesx.V <-  paste(namesx.v, rep(1:sv, each = length(namesx.v)), sep = "_")
    x.V <- XV[,which(namesxV %in% namesx.V), drop = F]
    }

    if(length(colnames.end.V) != 0) colnames.end.V <- paste(colnames.end.v,
                                                            rep(1:sv, each = length(colnames.end.v)), sep = "_")
    else colnames.end.V <- NULL
Hx.vne <- cbind(x.V, WxvD, instr.VV)

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
      nameswhV <- paste(nameswhv,rep(1:sv, each = ncol(whv)), sep = "_")
      hvD  <- ZV[, which(nameszv %in% nameswhV), drop = FALSE]
      WhvD <- matrix(0, ncol = ncol(hvD), nrow = n)
      seq_1 <- seq(1, ncol(hvD), sv)
      seq_2 <- seq(sv, ncol(hvD),  sv)
      for (i in 1: (ncol(hvD)/sv)) WhvD[,seq_1[i]:seq_2[i]] <-  as.matrix(Ws) %*% hvD[,seq_1[i]:seq_2[i]]
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



  endog <- list(colnames.end.f, colnames.end.V, col.end.f.l, colnames.end.vl)
  colinst <- c(nameInst, nameswhf, nameswhv)

  if(length(colinst) < length(c(colnames.end.f, colnames.end.V, col.end.f.l, colnames.end.vl)))
    stop("Not enough instruments specified: the model is not identified")
if((is.null(namesxfd) & is.null(namesxvd))) colinst  <- c("X", nameInst, nameswhf, nameswhv)
else   colinst <- c("X", "WX", nameInst, nameswhf, nameswhv)

}
  else{
    Hmat <- Zmat
    endog <- NULL
    colinst <-  NULL

  }



  ret <- list(y = y, Hmat = Hmat, Zmat = Zmat, endog = endog,
              instrum = colinst, l.split = l.split, Ws = Ws, nameswx, dur = dur)
  return(ret)
}

Matchnames <- function(nms, ct){
nms <- unlist(nms)
ct <- data.frame(id = letters[ct[,1]], gv = ct[,2])

place <- vector(mode = "list", length = nrow(ct))
hold <- vector(mode = "list", length = nrow(ct))
mhold <- vector(mode = "list", length = nrow(ct))
mholdr <- vector(mode = "list", length = nrow(ct))

for (i in 1: nrow(ct)) {
  place[[i]] <- which(endsWith(nms, paste0('_', match(ct$id[i] , letters) )))
  hold[[i]] <-  nms[place[[i]]]
  mhold[[i]] <- substr(hold[[i]], 1, nchar(hold[[i]])-2)
  mholdr[[i]] <- paste0(mhold[[i]], "_", ct$gv[i])
}

for(i in 1:nrow(ct)) nms[place[[i]]] <- mholdr[[i]]

  return(nms)
}

Matchgroups <- function(res, ct){

  ct <- data.frame(id = letters[ct[,1]], gv = ct[,2])

hold <- list()


place <- vector(mode = "list", length = nrow(ct))
hold <- vector(mode = "list", length = nrow(ct))
mhold <- vector(mode = "list", length = nrow(ct))
mholdr <- vector(mode = "list", length = nrow(ct))

for (i in 1: nrow(ct)) {
place[[i]] <- which(endsWith(rownames(res$coefficients), paste0('_', match(ct$id[i] , letters) )))
hold[[i]] <- rownames(res$coefficients)[place[[i]]]
mhold[[i]] <- substr(hold[[i]], 1, nchar(hold[[i]])-2)
mholdr[[i]] <- paste0(mhold[[i]], "_", ct$gv[i])
}

for(i in 1:nrow(ct)) rownames(res$coefficients)[place[[i]]] <- mholdr[[i]]


  return(res)
}
