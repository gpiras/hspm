
error_regimes <- function(formula, data, listw,  rgv, weps_rg = FALSE,
                          initial.value  = NULL, het,
                          verbose = FALSE, control, cl){




intro <- ols.data.prep.regimes(formula, data = data, listw = listw, rgv = rgv,
                               weps_rg = weps_rg, model = "error")


  y              <- as.matrix(intro[[1]])
  Hmat           <- intro[[2]]
  Zmat           <- intro[[3]]
  colnames.end   <- intro[[4]]
  colnames.instr <- intro[[5]]
  l.split        <- intro[[6]]
  Ws             <- intro[[7]]
  n              <- dim(Ws)[1]
  sv             <- l.split[[3]]
  nameswx        <- intro[[8]]
  dur            <- intro[[9]]
  ct             <- l.split[[6]]

  f.step <- spatial.ivreg.regimes(as.matrix(y), as.matrix(Zmat), as.matrix(Hmat), het)
  ubase  <- f.step[[3]]

### initial values for optimization
  pars <- in.val(weps_rg = weps_rg, initial.value = initial.value,
               Ws = Ws, ubase = ubase, sv = l.split[[3]])

##error part
  rhotilde <- error_part_regime(pars = pars, l.split = l.split,
                              Ws = Ws, het = het,
                              ubase = ubase, n = n, weps_rg = weps_rg,
                              verbose = verbose, control = control)

 out     <- co_transform(rhotilde = rhotilde, y = y,
                    Zmat = Zmat, Hmat = Hmat,
                    l.split = l.split, Ws = Ws, het = het)

 delta   <- out[[1]]
 utildeb <- out[[2]]

#print(delta)
#print(utildeb)
 res   <- error_efficient_regime(Ws = Ws, utildeb = utildeb,
                                n = n, weps_rg = weps_rg,
                                l.split = l.split, rhotilde = rhotilde,
                                Hmat = Hmat, Zmat = Zmat,
                                control = control, het = het,
                                verbose = verbose, delta = delta, y = y)

 res <- Matchgroups(res, ct)
 colnames.end <- Matchnames(colnames.end, ct)
 colnames.instr <- Matchnames(colnames.instr, ct)

 res    <- list(res, cl, colnames.end,  colnames.instr, nameswx, dur)
#print(res)
class(res) <- c("spregimes", "error_regimes")
return(res)
}

in.val <- function(weps_rg, initial.value = NULL, Ws, ubase, sv){
  if(weps_rg){
    if (is.null(initial.value)){
      Wubase <- Ws %*% ubase
      pars <- rep(coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1)), sv)
    }
    else {
      if(length(initial.value) != sv) pars <- rep(initial.value, sv)

    }
  }
  else{
    if (is.null(initial.value)){
      Wubase <- Ws %*% ubase
      pars <- coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1))
    }
    else {
      if(length(initial.value) != 1 )
        stop("Incorrect dimension of the the initial value")
      if(length(initial.value) != 1) pars <- initial.value
    }

  }
  return(pars)
}

error_part_regime <- function(pars, l.split,
                              Ws, het, ubase, n,
                              weps_rg, verbose, control){
  if(het){
    Ggmat <- gg_het_regime(Ws = Ws, u = ubase, n = n,
                           weps_rg = weps_rg, l.split = l.split)
    optres <- nlminb(pars, optimfunct_regime, lower= -0.9 + .Machine$double.eps ,
                     upper= 0.9 -  .Machine$double.eps, control= control,
                     v = Ggmat, weps_rg = weps_rg, verbose = verbose)
    rhotilde <- optres$par

  }
  else{

    Ggmat<-gg_hom_regime(Ws = Ws, u = ubase, n = n,
                         weps_rg = weps_rg, l.split = l.split)
    optres <- nlminb(pars, optimfunct_regime, lower = -0.9 + .Machine$double.eps ,
                     upper = 0.9 -  .Machine$double.eps, control = control,
                     v = Ggmat, weps_rg = weps_rg, verbose = verbose)
    rhotilde <- optres$par

  }
rhotilde
}

co_transform <- function(rhotilde, y, Zmat, Hmat, l.split, Ws, het){

  sv <- l.split[[3]]
  rgm <- l.split[[5]]

  if(length(rhotilde) == 1L){
    yt  <- y - rhotilde * Ws %*% y
    wZmat <- Ws %*% Zmat
    Zt <- Zmat - rhotilde * wZmat

    secondstep <- spatial.ivreg.regimes(y = yt , Zmat = Zt, Hmat = Hmat, het = het)
    delta <- coefficients(secondstep)
    utildeb <- y - Zmat %*% delta

  }
  else{

    yt  <- matrix(0, nrow = l.split[[1]], ncol = 1 )
    Zt    <- matrix(0, ncol = ncol(Zmat), nrow = nrow(Zmat))

    for(i in 1:sv) {

     yt[which(rgm[,i] == 1)]  <- ((y*rgm[,i]) - rhotilde[i] * Ws %*% (y*rgm[,i]))[which(rgm[,i] ==1)]
     Zt[,grep(paste("_", i, sep=""), colnames(Zmat))]   <- as.matrix(Zmat[,grep(paste("_", i, sep=""), colnames(Zmat))]) - rhotilde[i] * as.matrix(Ws %*% Zmat[,grep(paste("_", i, sep=""), colnames(Zmat))])

     }

    colnames(Zt) <- colnames(Zmat)

    secondstep <- spatial.ivreg.regimes(y = yt , Zmat = Zt, Hmat = Hmat, het = het)
    delta <- coefficients(secondstep)
    utildeb <- y - Zmat %*% delta

  }
  rownames(delta) <- colnames(Zmat)

  res <- list(delta = delta, utildeb = utildeb)
  return(res)
}

error_efficient_regime <- function(Ws, utildeb, n, weps_rg,
                                   l.split, rhotilde, Hmat,
                                   Zmat, control, het, verbose,
                                   delta, y){

  if(het){


    Ggmat <- gg_het_regime(Ws = Ws, u = utildeb, n = n,
                           weps_rg = weps_rg, l.split = l.split)
    gmm.weghts <- psirhorho_het_regime(rho = rhotilde,
                                       residuals = utildeb,
                                       Hmat = Hmat, Zmat = Zmat,
                                       Ws = Ws, l.split = l.split,
                                       weps_rg = weps_rg)

    optres <- nlminb(rhotilde, optimfunct_eff_regime, v = Ggmat,
                     vcmat = gmm.weghts$Phiinv, verbose = verbose,
                     lower = -0.9 + .Machine$double.eps ,
                     upper = 0.9 -  .Machine$double.eps,
                     control = control, weps_rg = weps_rg)

    rhofin <- optres$par
    gmm.weghts <- psirhorho_het_regime(rho = rhofin, residuals = utildeb,
                                       Hmat = Hmat, Zmat = Zmat,
                                       Ws = Ws, l.split = l.split,
                                       weps_rg = weps_rg)

    vcmat <- Omega_het_regime(rhofin, gmm.weghts$Pmat, gmm.weghts$A1,
                              gmm.weghts$A2, gmm.weghts$a.vec1,
                              gmm.weghts$a.vec2, Hmat,
                              Ggmat$bigG, gmm.weghts$Phiinv,
                              gmm.weghts$epsilon, gmm.weghts$Zstar, Ws,
                              l.split = l.split, weps_rg = weps_rg,
                              k = length(delta))




  }
  else{

    Ggmat <- gg_hom_regime(Ws = Ws, u = utildeb, n = n,
                           weps_rg = weps_rg, l.split = l.split)

    gmm.weghts <- psirhorho_hom_regime(rho = rhotilde, residuals = utildeb,
                                       Hmat = Hmat, Zmat = Zmat, Ws = Ws,
                                       d = Ggmat$d, v.vec = Ggmat$v.vec,
                                       l.split = l.split, weps_rg = weps_rg)
    #print(utildeb)
    optres <- nlminb(rhotilde, optimfunct_eff_regime, v = Ggmat,
                     vcmat = gmm.weghts$Phiinv, verbose = verbose,
                     control = control, lower = -1 + .Machine$double.eps ,
                     upper = 1 -  .Machine$double.eps,
                     weps_rg = weps_rg)

    rhofin <- optres$par
    #print(rhofin)
    gmm.weghts <-psirhorho_hom_regime(rho = rhofin, residuals = utildeb,
                                      Hmat = Hmat, Zmat = Zmat, Ws = Ws,
                                      d = Ggmat$d, v.vec = Ggmat$v.vec,
                                      l.split = l.split, weps_rg = weps_rg)
    #   print(utildeb)
    vcmat <- Omega_hom_regime(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2,
                              gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat,
                              Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar,
                              l.split, weps_rg, k = length(delta))




  }

  if(weps_rg){
    k <- length(delta)
    sv <- l.split[[3]]
    svm <- l.split[[4]]
    seq1 <- seq(1, k+sv, k/sv+1)
    seq2 <- seq(k/sv+1, k+sv, k/sv+1)
    seq3 <- seq(1, k+sv, k/sv)
    seq4 <- seq(k/sv, k+sv, k/sv)

    deltarho1 <- vector(mode = "numeric", length = sv+k)

    for(i in 1:sv) {
      deltarho1[seq1[i]:seq2[i]] <- c(delta[grep(paste("_", i, sep=""), rownames(delta))], rhotilde[i])
      names(deltarho1)[seq1[i]:seq2[i]] <- c(rownames(delta)[grep(paste("_", i, sep=""), rownames(delta))],paste("We_", i, sep = "") )
    }

  }
  else{
    deltarho1 <- c(as.numeric(delta), rhofin)
    names(deltarho1) <- c(rownames(delta), "We")
  }

  deltarho <- matrix(deltarho1, ncol = 1, nrow = length(deltarho1))
  rownames(deltarho) <- names(deltarho1)
   yp <- Zmat %*% delta
   yp <- array(yp, dim = c(length(yp), 1),
              dimnames = list(seq(1, length(yp)), ""))
   utildeb <- array(utildeb, dim = c(length(utildeb), 1),
               dimnames = list(seq(1, length(utildeb)), ""))

   res <- list(coefficients = deltarho, var = vcmat,
                residuals = utildeb, X = Zmat,
                y = y, yp = yp)

  return(res)
}
