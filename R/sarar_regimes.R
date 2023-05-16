sarar_regimes <- function(formula, data, listw,  rgv, het,
                          weps_rg = weps_rg, wy_rg = wy_rg,
                          initial.value  = NULL,  verbose = FALSE, control, cl){


  intro <- iv.lag.data.prep.regimes(formula, data = data, listw = listw,
                                    wy_rg = wy_rg, rgv = rgv, weps_rg = weps_rg)

  y        <- as.matrix(intro[[1]])
  Hmat     <- intro[[2]]
  Zmat     <- intro[[3]]
  colnames.end   <- intro[[4]]
  colnames.instr <- intro[[5]]
  l.split <- intro[[6]]
  Ws <- intro[[7]]
  colinstr <- intro[[8]]
  n <- dim(Ws)[1]
  sv <- l.split[[3]]
  ct       <- intro[[6]][[6]]


  f.step <- spatial.ivreg.regimes(as.matrix(y), as.matrix(Zmat), as.matrix(Hmat), het)
  ubase <- f.step[[3]]

  ### initial values for optimization
  pars <- in.val(weps_rg = weps_rg, initial.value = initial.value,
                 Ws = Ws, ubase = ubase, sv = sv)

  ##error part
  rhotilde <- error_part_regime(pars = pars, l.split = l.split,
                                Ws = Ws, het = het,
                                ubase = ubase, n = n, weps_rg = weps_rg,
                                verbose = verbose, control = control)

  ##co_transform_sarar includes second estimation
  out <- co_transform_sarar(rhotilde = rhotilde, y = y,
                      Zmat = Zmat, Hmat = Hmat,
                      l.split = l.split, Ws = Ws, het = het, wy_rg = wy_rg)

  delta <- out[[1]]
  utildeb <- out[[2]]
  ##calculates the vc matrix final

  res <- error_efficient_regime(Ws = Ws, utildeb = utildeb,
                                n = n, weps_rg = weps_rg,
                                l.split = l.split,
                                rhotilde = rhotilde,
                                Hmat = Hmat, Zmat = Zmat,
                                control = control, het = het,
                                verbose = verbose, delta = delta, y = y)

  res <- Matchgroups(res, ct)
  colnames.end <- Matchnames(colnames.end, ct)
  colnames.instr <- Matchnames(colnames.instr, ct)

  res <- list(res, cl, colnames.end,  colnames.instr, colinstr)
  #print(res)
  class(res) <- c("spregimes", "sarar_regimes")
  return(res)
}



in.val <- function(weps_rg, initial.value = NULL, Ws, ubase, sv){
  if(weps_rg){
    if (is.null(initial.value)){
      Wubase <- Ws %*% as.matrix(ubase)
      pars <- rep(coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1)), sv)
    }
    else {
      if(length(initial.value) != sv) pars <- rep(initial.value, sv)

    }
  }
  else{
    if (is.null(initial.value)){
      Wubase <- Ws %*% as.matrix(ubase)
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

# error_part_regime <- function(pars, l.split,
#                               Ws, het, ubase, n,
#                               weps_rg, verbose, control){
#   if(het){
#     Ggmat <- gg_het_regime(Ws = Ws, u = ubase, n = n,
#                            weps_rg = weps_rg, l.split = l.split)
#     optres <- nlminb(pars, optimfunct_regime, lower= -0.9 + .Machine$double.eps ,
#                      upper= 0.9 -  .Machine$double.eps, control= control,
#                      v = Ggmat, weps_rg = weps_rg, verbose = verbose)
#     rhotilde <- optres$par
#
#   }
#   else{
#
#     Ggmat<-gg_hom_regime(Ws = Ws, u = ubase, n = n,
#                          weps_rg = weps_rg, l.split = l.split)
#     optres <- nlminb(pars, optimfunct_regime, lower = -0.9 + .Machine$double.eps ,
#                      upper = 0.9 -  .Machine$double.eps, control = control,
#                      v = Ggmat, weps_rg = weps_rg, verbose = verbose)
#     rhotilde <- optres$par
#
#   }
#   rhotilde
# }
#
co_transform_sarar <- function(rhotilde, y, Zmat, Hmat, l.split, Ws, het, wy_rg){

  if(length(rhotilde) == 1L){
    yt  <- y - rhotilde * Ws %*% y
    wZmat <- Ws %*% Zmat
    Zt <- Zmat - rhotilde * wZmat

    secondstep <- spatial.ivreg.regimes(y = yt , Zmat = Zt, Hmat = Hmat, het = het)
    delta <- coefficients(secondstep)
    utildeb <- y - Zmat %*% delta
  }
  else{
    if(!wy_rg) stop("if weps_rg is TRUE also wy_rg should be TRUE")
    sv <- l.split[[3]]
    rgm <- l.split[[5]]

    yt  <- matrix(0, nrow = l.split[[1]], ncol = 1 )
    for(i in 1:sv) yt[which(rgm[,i] == 1)]  <- ((y*rgm[,i]) - rhotilde[i] * Ws %*% (y*rgm[,i]))[which(rgm[,i] ==1)]
    #this multiplies each colums of zmat for the corresponding rho
    Zt    <- matrix(0, ncol = ncol(Zmat), nrow = nrow(Zmat))
    for(i in 1:sv) Zt[,grep(paste("_", i, sep=""), colnames(Zmat))]   <- as.matrix(Zmat[,grep(paste("_", i, sep=""), colnames(Zmat))]) - rhotilde[i] * as.matrix(Ws %*% Zmat[,grep(paste("_", i, sep=""), colnames(Zmat))])
    colnames(Zt) <- colnames(Zmat)

     secondstep <- spatial.ivreg.regimes(y = yt , Zmat = Zt, Hmat = Hmat, het = het)
    delta <- coefficients(secondstep)
    utildeb <- y - Zmat %*% delta
  }
  rownames(delta) <- colnames(Zmat)
 #print(delta)
  out <- list(delta = delta, utildeb = utildeb)
  return(out)
}
