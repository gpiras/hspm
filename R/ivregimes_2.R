

##### Functions for regimes ####
#' @title Estimation of spatial regime models with endogenous variables
#' @name ivregimes
#' @param formula a symbolic description of the model of the form \code{y ~ x_f | x_v | h_f | h_v} where \code{y} is the dependent variable, \code{x_f} are the regressors that do not vary by regimes,  \code{x_v} are the regressors that vary by regimes, \code{h_f} are the fixed instruments and \code{h_v} are the instruments that vary by regimes.
#' @param data the data of class \code{data.frame}.
#' @param rgv an object of class \code{formula} to identify the regime variables
#' @param vc   one of \code{c('classical', 'robust', 'OGMM')}. If \code{OGMM} an optimal weighted GMM is used to estimate the VC matrix (for additional details see Anselin and Rey, 2014).
#' @param object an object of class ivregime
#' @param ... additional arguments
#' @param x an object of class ivregime
#' @param digits number of digits
#'
#'
#' @details
#'
#' The model estimated is:
#'
#' \deqn{
#' y_{ij}= \mathbf{x_{ij,k}}\beta_j + \mathbf{Y_{ij,k}}\gamma_j + \epsilon
#' }
#' for i=1,..,n representing the sample observations, and j = 1,..., J representing
#' the  regimes

#' @examples
#' data("natreg")
#' form   <- HR90  ~ 0 | MA90 + PS90 + RD90 + UE90 | 0 | MA90 + PS90 + RD90 + FH90 + FP89 + GI89
#' split  <- ~ REGIONS
#' mod <- ivregimes(formula = form, data = natreg, rgv = split, vc = "robust")
#' summary(mod)
#' mod1 <- ivregimes(formula = form, data = natreg, rgv = split, vc = "OGMM")
#' summary(mod1)
#' form1   <- HR90  ~ MA90 + PS90 |  RD90 + UE90 -1 | MA90 + PS90 | RD90 + FH90 + FP89 + GI89 -1
#' mod2 <- ivregimes(formula = form1, data = natreg, rgv = split, vc = "classical")
#' summary(mod2)
#'
#' @author Gianfranco Piras and Mauricio Sarrias
#' @return An object of class \code{ivregimes}
#' @import Formula sphet stats
#' @export

ivregimes <- function(formula, data, rgv = NULL,
                      vc = c("classical", "robust", "OGMM")){

  cl <- match.call()



  #Obtain arguments
  vc                <- match.arg(vc)

  #process the data
  intro    <- iv.data.prep.regimes(formula = formula, data = data,
                                   rgv = rgv)

  #split the object list intro:
  # data
  # formula
  # k tot number of variables in model
  # totr tot of variables varying by regime
  y        <- intro[[1]]
  Hmat     <- intro[[2]]
  Zmat     <- intro[[3]]
  colnames.end.f <- intro[[4]]
  colnames.end.v <- intro[[5]]
  colnames.instr <- intro[[6]]

  res <- tsls_regimes(y, Hmat, Zmat, vc)

  res <- list(res, cl, colnames.end.f, colnames.end.v, colnames.instr)
  class(res) <- "ivregimes"
  return(res)
}


tsls_regimes <- function(y, Hmat, Zmat, vc){


  df <- nrow(Zmat) - ncol(Zmat)
  n <- nrow(Zmat)
  HH <- crossprod(Hmat)
  Hye <- crossprod(Hmat, Zmat)
  bz <- solve(HH, Hye)
  Zp <- Hmat %*% bz
  ZpZpi <- solve(crossprod(Zp))
  betaiv <- ZpZpi %*% crossprod(Zp, as.matrix(y))

  yp <- Zmat %*% betaiv
  e  <- as.matrix(y) - yp


  ##simple
  if(vc =="classical"){

    cpe <- crossprod(e)
    vcmatrix <- (as.numeric(cpe) /df) * solve(crossprod(Zp))

    return(list(betaiv, vcmatrix))
  }


  if(vc == "robust"){
    e2 <- e^2
    ZoZ <- crossprod(Zp, (as.matrix(Zp) * as.numeric(e2)))
    vcmatrix <-  ZpZpi %*% ZoZ %*% ZpZpi * (n/df)
    return(list(betaiv, vcmatrix))

  }


  if(vc == "OGMM"){
    e2 <- e^2
    uH <- matrix(0, nrow = dim(Hmat)[1], ncol = dim(Hmat)[2])
    for (i in 1 : dim(Hmat)[2]) uH[,i] <- e2 * Hmat[,i]
    Smat <- crossprod(Hmat, uH)
    Smati <- solve(Smat)
    ZpH <- crossprod(Zmat, Hmat)
    HpZ <- crossprod(Hmat, Zmat)
    fp <- solve(ZpH %*% Smati %*% HpZ)
    Hpy <- crossprod(Hmat, as.matrix(y))
    sp <- ZpH %*% Smati %*% Hpy
    b_OWGMM <- fp %*% sp
    vcmatrix <-  solve(ZpH %*% Smati %*% HpZ)

    return(list(b_OWGMM, vcmatrix))
  }

}


iv.data.prep.regimes <- function(formula, data, rgv){
  n                <- dim(data)[1]
  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- as.numeric(unique(splitvar))
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))

  F1 <- Formula(formula)

  parts <- length(F1)

  if(parts[2] != 4) stop("Formula should have four parts")

  mf <- model.frame(F1, data = data)

  y <- model.part(F1, data = mf, lhs = 1, drop = FALSE)


  ## extract X fixed
  Xf <- model.matrix(F1, data = mf, rhs = 1, drop = FALSE)
  namesxf <- colnames(Xf)
  ## extract X variable
  Xv <- model.matrix(F1, data = mf, rhs = 2, drop = FALSE)
  namesxv <- colnames(Xv)
  namesZ <- c(namesxf, namesxv)

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
  namesxv <- paste(namesxv,rep(1:sv,each = k2), sep = "_")
  colnames(XV) <- namesxv

  ## extract instrument
  Zf <- model.matrix(F1, data = mf, rhs = 3, drop = FALSE)
  nameszf <- colnames(Zf)
  Zv <- model.matrix(F1, data = mf, rhs = 4, drop = FALSE)
  nameszv <- colnames(Zv)
  namesH <- c(nameszf, nameszv)
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

  ## if instruments check identification
  end.f <- Xf[, !(colnames(Xf) %in% colnames(Zf)), drop = FALSE]
  end.v <- Xv[, !(colnames(Xv) %in% colnames(Zv)), drop = FALSE]
  colnames.end.f <- colnames(end.f)
  colnames.end.v <- colnames(end.v)
  # cat("the endogenous fixed variables are ", colnames(end.f), "\n")
  # cat("and the instruments are ", colnames(Zf), "\n")
  #
  # cat("the endogenous split variables are ", colnames(end.v), "\n")
  # cat("and the instruments are ", colnames(Zv), "\n")


  Hmat <- cbind(as.matrix(Zf), as.matrix(ZV))
  Zmat <- cbind(as.matrix(Xf), as.matrix(XV))
  colinst <- namesH[-which((namesH %in% namesZ))]

  if(length(colinst) < length(c(colnames.end.f, colnames.end.v)))
    stop("Not enough instruments specified: the model is not identified")

  ret <- list(y = y, Hmat = Hmat, Zmat = Zmat, colnames.end.f, colnames.end.v, colinst)
  return(ret)
}



### S3 methods ----

#' @rdname ivregimes
#' @method coef ivregimes
#' @export
coef.ivregimes <- function(object, ...){
  object[[1]][[1]]
}


#' @rdname ivregimes
#' @method vcov ivregimes
#' @import stats
#' @export
vcov.ivregimes <- function(object, ...){
  V <- object[[1]][[2]]
  return(V)
}


#' @rdname ivregimes
#' @method print ivregimes
#' @import stats
#' @export
print.ivregimes <- function(x,
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



#' @rdname ivregimes
#' @method summary ivregimes
#' @import stats
#' @export
summary.ivregimes <- function(object, ...){
  b                   <- coef(object)
  std.err             <- sqrt(diag(vcov(object)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.ivregimes", "ivregimes")
  return(object)
}


#' @rdname ivregimes
#' @method print summary.ivregimes
#' @import stats
#' @export
print.summary.ivregimes <- function(x,
                                    digits = max(5, getOption("digits") - 3),
                                    ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                          IV Regimes Model \n")
  cat("        ------------------------------------------------------------\n")
   cat("\nCall:\n")
   cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

   cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

  cat("\nEndogenous variables:\n")

  cat(paste(x[[3]], x[[4]], sep=" "))

  cat("\nInstruments:\n")

  cat(paste(x[[5]], sep=" "))
  invisible(x)
}

