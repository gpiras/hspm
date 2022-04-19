##### Functions for linear model with lagged explanatory variables####
#' Estimation of spatial regime models
#' @name error_regimes
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{data.frame}.
#' @param listw a spatial weighting matrix
#' @param initial.value startingo point for the optimization
#' @param rgv variable to identify the regimes
#' @param het heteroskedastic variance-covariance matrix
#' @param cl record calls
#' @param object an object of class lag_regime
#' @param ... additional arguments
#' @param x an object of class lag_regime
#' @param digits number of digits

#'
error_regimes <- function(formula, data, listw,  rgv, weps_rg,
                          initial.value, het,
                          step1.c, control, cl){




intro <- ols.data.prep.regimes(formula, data = data, listw = listw, rgv = rgv)


  y        <- intro[[1]]
  Hmat     <- intro[[2]]
  Zmat     <- intro[[3]]
  colnames.end<- intro[[4]]
  colnames.instr <- intro[[5]]
  l.split <- intro[[6]]
  Ws <- intro[[7]]


f.step <- spatial.ivreg.regimes(as.matrix(y), as.matrix(Zmat), as.matrix(Hmat), het)
ubase <- f.step[[3]]

### initial values for optimization
pars <- in.val(weps_rg, initial.value, Ws, ubase)




  return(res)

}

### S3 methods ----

#' @rdname error_regimes
#' @method coef error_regimes
#' @export
coef.error_regimes <- function(object, ...){
  object[[1]][[1]]
}


#' @rdname error_regimes
#' @method vcov error_regimes
#' @import stats
#' @export
vcov.ols_regimes <- function(object, ...){
  V <- object[[1]][[2]]
  return(V)
}




#' @rdname error_regimes
#' @method print error_regimes
#' @import stats
#' @export
print.error_regimes <- function(x,
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



#' @rdname error_regimes
#' @method summary error_regimes
#' @import stats
#' @export
summary.error_regimes <- function(object, ...){
  b                   <- coef(object)
  std.err             <- sqrt(diag(vcov(object)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.error_regimes", "error_regimes")
  return(object)
}


#' @rdname error_regimes
#' @method print summary.error_regimes
#' @import stats
#' @export
print.summary.error_regimes <- function(x,
                                      digits = max(5, getOption("digits") - 3),
                                      ...)
{
  if(!is.null(x[[4]])){
    cat("        ------------------------------------------------------------\n")
    cat("                         Spatial Error Regimes Model      \n")
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
    cat("                       Spatial Error Regimes Model       \n")
    cat("        ------------------------------------------------------------\n")
    cat("\nCall:\n")
    cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    cat("\nCoefficients:\n")
    printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)


  }


  invisible(x)
}







in.val <- function(weps_rg, initial.value, Ws, ubase){
  if(weps_rg){
    if (is.null(initial.value)){
      Wubase <- Ws %*% ubase
      pars <- rep(coefficients(lm(as.numeric(ubase) ~ as.numeric(Wubase)-1)),2)
    }
    else {
      if(length(initial.value) != 1 || length(initial.value) != 2)
        stop("Incorrect dimension of the the initial values")
      if(length(initial.value) != 2) pars <- rep(initial.value, 2)
      else pars <- initial.values
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
      if(length(initial.value) != 1) pars <- initial.values
    }

  }
  return(pars)
}


error_part_regime <- function(pars, l.split, het, ubase, n, weps_rg){
  if(het){
    Ggmat <- gg_het_regime(Ws, ubase, n, weps_rg, l.split)


    optres <- nlminb(pars, optimfunct_regime, lower= -0.9 + .Machine$double.eps ,
                     upper= 0.9 -  .Machine$double.eps, control= control,
                     v = Ggmat)
    rhotilde<-optres$par


    if(step1.c){
      gmm.weghts1.c <- psirhorho_het_regime(rhotilde, ubase, Hmat, Zmat, Ws, step1.c = TRUE)

      optres <- nlminb(rhotilde, optimfunct_eff_regime, v = Ggmat, vcmat = gmm.weghts1.c$Phiinv,
                       verbose = verbose, lower = -0.9 + .Machine$double.eps,
                       upper = 0.9 -  .Machine$double.eps, control = control)
      rhotilde <- optres$par
      gmm.weghts1.c <- psirhorho_het_regime(rhotilde, ubase, Hmat, Zmat, Ws, step1.c = TRUE)
      vcmat_2sls <- Omega_het_regime(rhotilde, gmm.weghts1.c$Pmat, gmm.weghts1.c$A1,
                              gmm.weghts1.c$A2, gmm.weghts1.c$a.vec1,
                              gmm.weghts1.c$a.vec2, Hmat, Ggmat$bigG,
                              gmm.weghts1.c$Phiinv, gmm.weghts1.c$epsilon,
                              gmm.weghts1.c$Zstar, Ws, step1.c = TRUE)


      coeff_2sls <- as.matrix(c(coefficients(firststep), rhotilde))
      rownames(coeff_2sls)<-c(colnames(Zmat), 'Wu')
      s2_2sls<-crossprod(ubase)/(n-k)



      results_2sls <- list(coefficients = coeff_2sls, var = vcmat_2sls$Omega,
                           residuals = as.numeric(ubase),
                           firststep = firststep$coefficients, init.rho = rhotilde)

    }


  }
  else{

    Ggmat<-gg_hom(Ws, ubase, n)
    optres <- nlminb(pars, optimfunct, control = control,
                     v = Ggmat, verbose = verbose, lower = -0.9 + .Machine$double.eps ,
                     upper = 0.9 -  .Machine$double.eps)
    rhotilde <- optres$par

  }

}
