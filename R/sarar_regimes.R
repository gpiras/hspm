##### Functions for linear model with spatially lagged dependent variable, error dependence and regimes####
#' Estimation of spatial regime models
#' @name sarar_regimes
#' @param formula a symbolic description of the model of the form \code{y ~ x_f | x_v | wx | h_f | h_v | wh} where \code{y} is the dependent variable, \code{x_f} are the regressors that do not vary by regimes,  \code{x_v} are the regressors that vary by regimes, \code{wx} are the spatially lagged regressors, \code{h_f} are the instruments that do not vary by regimes,  \code{h_v} are the instruments that vary by regimes, \code{wh} are the spatially lagged instruments.
#' @param data the data of class \code{data.frame}.
#' @param listw a spatial weighting matrix of class \code{listw}, \code{matrix} or \code{Matrix}
#' @param wy_rg default \code{wy_rg = FALSE}, the lagged dependent variable does not vary by regime (see details)
#' @param weps_rg default FALSE, if TRUE the spatial error term varies by regimes (see details)
#' @param initial.value initial value for the spatial error parameter
#' @param rgv an object of class \code{formula} to identify the regime variables
#' @param het heteroskedastic variance-covariance matrix
#' @param cl record calls
#' @param verbose print a trace of the optimization
#' @param control argument for optimization
#' @param object an object of class sarar_regimes
#' @param ... additional arguments
#' @param x an object of class sarar_regimes
#' @param digits number of digits
#'
#'
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
  Ws <- listw
  n <- dim(Ws)[1]
  sv <- l.split[[3]]

  f.step <- spatial.ivreg.regimes(as.matrix(y), as.matrix(Zmat), as.matrix(Hmat), het)
  ubase <- f.step[[3]]

  ### initial values for optimization
  pars <- in.val(weps_rg = weps_rg, initial.value = initial.value,
                 Ws = Ws, ubase = ubase, sv = l.split[[3]])

  ##error part
  rhotilde <- error_part_regime(pars = pars, l.split = l.split,
                                Ws = Ws, het = het,
                                ubase = ubase, n = n, weps_rg = weps_rg,
                                verbose = verbose, control = control)

  out <- co_transform_sarar(rhotilde = rhotilde, y = y,
                      Zmat = Zmat, Hmat = Hmat,
                      l.split = l.split, Ws = Ws, het = het, wy_rg = wy_rg)

  delta <- out[[1]]
  utildeb <- out[[2]]
  #k <- out[[3]]
  #print(delta)
  res <- error_efficient_regime(Ws = Ws, utildeb = utildeb,
                                n = n, weps_rg = weps_rg,
                                l.split = l.split,
                                rhotilde = rhotilde,
                                Hmat = Hmat, Zmat = Zmat,
                                control = control, het = het,
                                verbose = verbose, delta = delta)

  res <- list(res, cl, colnames.end,  colnames.instr)
  #print(res)
  class(res) <- "sarar_regimes"
  return(res)
}

### S3 methods ----

#' @rdname sarar_regimes
#' @method coef sarar_regimes
#' @export
coef.sarar_regimes <- function(object, ...){
  object[[1]][[1]]
}


#' @rdname sarar_regimes
#' @method vcov sarar_regimes
#' @import stats
#' @export
vcov.sarar_regimes <- function(object, ...){
  V <- object[[1]][[2]]
  return(V)
}




#' @rdname sarar_regimes
#' @method print sarar_regimes
#' @import stats
#' @export
print.sarar_regimes <- function(x,
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



#' @rdname sarar_regimes
#' @method summary sarar_regimes
#' @import stats
#' @export
summary.sarar_regimes <- function(object, ...){
  b                   <- coef(object)
  std.err             <- sqrt(diag(vcov(object)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.sarar_regimes", "sarar_regimes")
  return(object)
}


#' @rdname sarar_regimes
#' @method print summary.sarar_regimes
#' @import stats
#' @export
print.summary.sarar_regimes <- function(x,
                                        digits = max(5, getOption("digits") - 3),
                                        ...)
{
  if(!is.null(x[[4]])){
    cat("        ------------------------------------------------------------\n")
    cat("                         Spatial SARAR Regimes Model      \n")
    cat("                      and additional endogenous variables               \n")
    cat("        ------------------------------------------------------------\n")
    cat("\nCall:\n")
    cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    cat("\nCoefficients:\n")
    printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)


    cat("\nEndogenous variables:\n")

    cat(paste(unlist(x[[3]]), sep=" "))

    cat("\nInstruments:\n")

    cat(paste(x[[4]], sep=" "))
  }

  else{
    cat("        ------------------------------------------------------------\n")
    cat("                       Spatial SARAR Regimes Model       \n")
    cat("        ------------------------------------------------------------\n")
    cat("\nCall:\n")
    cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    cat("\nCoefficients:\n")
    printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)


  }


  invisible(x)
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

error_part_regime <- function(pars, l.split, Ws, het, ubase, n, weps_rg, verbose, control){
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
    for(i in 1: sv) Zt[,grep(paste("_", i, sep=""), colnames(Zmat))]   <- as.matrix(Zmat[,grep(paste("_", i, sep=""), colnames(Zmat))]) - rhotilde[i] * as.matrix(Ws %*% Zmat[,grep(paste("_", i, sep=""), colnames(Zmat))])
    colnames(Zt) <- colnames(Zmat)
    secondstep <- spatial.ivreg.regimes(y = yt , Zmat = Zt, Hmat = Hmat, het = het)
    delta <- coefficients(secondstep)
    utildeb <- y - Zmat %*% delta
  }
  rownames(delta) <- colnames(Zmat)

  out <- list(delta = delta, utildeb = utildeb)
  return(out)
}

#
# #error_efficient_regime <- function(Ws, utildeb, n, weps_rg, l.split,
#  #                                  rhotilde, Hmat,
#   #                                 Zmat, control, het,
#    #                                verbose, delta){
#   if(het){
#
#
#     Ggmat <- gg_het_regime(Ws = Ws, u = utildeb, n = n,
#                            weps_rg = weps_rg, l.split = l.split)
#     gmm.weghts <- psirhorho_het_regime(rho = rhotilde,
#                                        residuals = utildeb,
#                                        Hmat = Hmat, Zmat = Zmat,
#                                        Ws = Ws, l.split = l.split,
#                                        weps_rg = weps_rg)
#
#     optres <- nlminb(rhotilde, optimfunct_eff_regime, v = Ggmat,
#                      vcmat = gmm.weghts$Phiinv, verbose = verbose,
#                      lower = -0.9 + .Machine$double.eps ,
#                      upper = 0.9 -  .Machine$double.eps,
#                      control = control, weps_rg = weps_rg)
#
#     rhofin <- optres$par
#     gmm.weghts <- psirhorho_het_regime(rho = rhofin, residuals = utildeb,
#                                        Hmat = Hmat, Zmat = Zmat,
#                                        Ws = Ws, l.split = l.split,
#                                        weps_rg = weps_rg)
#
#     vcmat <- Omega_het_regime(rhofin, gmm.weghts$Pmat, gmm.weghts$A1,
#                               gmm.weghts$A2, gmm.weghts$a.vec1,
#                               gmm.weghts$a.vec2, Hmat,
#                               Ggmat$bigG, gmm.weghts$Phiinv,
#                               gmm.weghts$epsilon, gmm.weghts$Zstar, Ws,
#                               l.split = l.split, weps_rg = weps_rg,
#                               k = length(delta))
#
#
#
#
#   }
#   else{
#
#     Ggmat <- gg_hom_regime(Ws = Ws, u = utildeb, n = n,
#                            weps_rg = weps_rg, l.split = l.split)
#
#     gmm.weghts <- psirhorho_hom_regime(rho = rhotilde, residuals = utildeb,
#                                        Hmat = Hmat, Zmat = Zmat, Ws = Ws,
#                                        d = Ggmat$d, v.vec = Ggmat$v.vec,
#                                        l.split = l.split, weps_rg = weps_rg)
#     #print(utildeb)
#     optres <- nlminb(rhotilde, optimfunct_eff_regime, v = Ggmat,
#                      vcmat = gmm.weghts$Phiinv, verbose = verbose,
#                      control = control, lower = -1 + .Machine$double.eps ,
#                      upper = 1 -  .Machine$double.eps,
#                      weps_rg = weps_rg)
#
#     rhofin <- optres$par
#     #print(utildeb)
#     gmm.weghts <-psirhorho_hom_regime(rho = rhofin, residuals = utildeb,
#                                       Hmat = Hmat, Zmat = Zmat, Ws = Ws,
#                                       d = Ggmat$d, v.vec = Ggmat$v.vec,
#                                       l.split = l.split, weps_rg = weps_rg)
#     #   print(utildeb)
#     vcmat <- Omega_hom_regime(rhofin, gmm.weghts$Pmat, gmm.weghts$A1, gmm.weghts$A2,
#                               gmm.weghts$a.vec1, gmm.weghts$a.vec2, Hmat,
#                               Ggmat$bigG, gmm.weghts$Phiinv, gmm.weghts$epsilon, gmm.weghts$Zstar,
#                               l.split, weps_rg, k = length(delta))
#
#
#
#
#   }
#
#   if(weps_rg){
#     k <- length(delta)
#     sv <- l.split[[3]]
#     svm <- l.split[[4]]
#     seq1 <- seq(1, k+sv, k/sv+1)
#     seq2 <- seq(k/sv+1, k+sv, k/sv+1)
#     seq3 <- seq(1, k+sv, k/sv)
#     seq4 <- seq(k/sv, k+sv, k/sv)
#
#     deltarho <- vector(mode = "numeric", length = sv+k)
#
#     for(i in 1:sv) {
#       deltarho[seq1[i]:seq2[i]] <- c(delta[grep(paste("_", i, sep=""), rownames(delta))], rhotilde[i])
#       names(deltarho)[seq1[i]:seq2[i]] <- c(rownames(delta)[grep(paste("_", i, sep=""), rownames(delta))],paste("We_", i, sep = "") )
#     }
#
#   }
#   else{
#     deltarho <- c(delta, rhofin)
#     names(deltarho) <- c(rownames(delta), "We")
#   }
#   result <- list(coefficients = deltarho, var = vcmat, residuals = utildeb)
#
#   return(result)
#
# }
# #
