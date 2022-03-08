


##### Functions for regimes ####
#' Estimation of spatial regime models
#' @name ivregimes
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{data.frame}.
#' @param rgv variable to identify the regimes
#' @param vc   one of ("classical", "robust", "OGMM")
#' @param endog formula containing  endogenous variables
#' @param endog_v logical vector defining which endogenous variable should vary by regime
#' @param instr formula containing instrument for endogenous variables
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
#' for i=1,..,n representing the sample observations, and j =1,..., J representing
#' the  regimes
#'
#'
#' @author Gianfranco Piras and Mauricio Sarrias
#' @return An object of class ``\code{spregimes}''
#' @import Formula sphet stats
#' @export

ivregimes <- function(formula, data, rgv = NULL,
                      endog = NULL, instr = NULL,
                      endog_v = NULL,
                    vc = c("classical", "robust", "OGMM")){

  cl <- match.call()

  if(is.null(endog)) stop("Endogenous variables not specified, use function regimes")
  if(is.null(instr)) stop("Instruments not specified for endogenous variables")

  if(is.null(rgv)) stop("regimes variable not specified")
  if(class(rgv) != "formula") stop("regimes variable has to be a formula")

  if(is.null(endog_v)) endog_v <- rep(TRUE, length(all.vars(endog)))

   #Obtain arguments
  vc                <- match.arg(vc)

  #process the data
  intro    <- iv.data.prep.regimes(formula = formula, data = data,
                                   rgv = rgv, endog = endog,
                                   instr = instr, endog_v = endog_v)

  #split the object list intro:
  # data
  # formula
  # k tot number of variables in model
  # totr tot of variables varying by regime
  y        <- intro[[1]]
  Hmat     <- intro[[2]]
  Zmat     <- intro[[3]]

  res <- tsls_regimes(y, Hmat, Zmat, vc)

 # out <- list(res, cl)
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
  betaiv <- ZpZpi %*% crossprod(Zp, y)

  yp <- Zmat %*% betaiv
  e <- y - yp


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
    uH <- matrix(0, nrow = dim(Hmat)[1], ncol = dim(Hmat)[2])
    for (i in 1 : dim(Hmat)[2]) uH[,i] <- e2 * Hmat[,i]
    Smat <- crossprod(Hmat, uH)
    Smati <- solve(Smat)
    ZpH <- crossprod(Zmat, Hmat)
    HpZ <- crossprod(Hmat, Zmat)
    fp <- solve(ZpH %*% Smati %*% HpZ)
    Hpy <- crossprod(Hmat, y)
    sp <- ZpH %*% Smati %*% Hpy
    b_OWGMM <- fp %*% sp
    vcmatrix <-  solve(ZpH %*% Smati %*% HpZ)

      return(list(b_OWGMM, vcmatrix))
  }

}


iv.data.prep.regimes <- function(formula, data, rgv, endog,
                                 instr, endog_v){
  n                <- dim(data)[1]
  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- as.numeric(unique(splitvar))
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))

  endog            <- as.matrix(lm(endog, data, method="model.frame"))
  if (length(endog_v) != ncol(endog))
    stop("the number of endogenous variables differ from endog_v")

  endogv         <- as.matrix(endog[,which(endog_v == TRUE)])
  endogf         <- as.matrix(endog[,which(endog_v == FALSE)])
  namesev        <- colnames(endog)[which(endog_v == TRUE)]
  namesef        <- colnames(endog)[which(endog_v == FALSE)]
  instr          <- as.matrix(lm(instr, data, method="model.frame"))
  namesinst      <- colnames(instr)

  #### all instruments vary
  p              <- dim(instr)[2]
  instr_r        <- matrix(0, ncol = 2*p, nrow = n)
  seq_1          <- seq(1, sv*p, p)
  seq_2          <- seq(p, sv*p,  p)
  for(i in 1:sv) instr_r[,seq_1[i]:seq_2[i]] <- instr * rgm[,i]
  namesinstr_r <- paste(namesinst,rep(1:sv, each = p), sep = "_")

  ### only endogv vary
  q              <- dim(endogv)[2]
  endogvv        <- matrix(0, ncol = 2*q, nrow = n)
  seq_1          <- seq(1, sv*q, q)
  seq_2          <- seq(q, sv*q,  q)
  for(i in 1:sv) endogvv[,seq_1[i]:seq_2[i]] <- endogv * rgm[,i]
  namesev_r     <- paste(namesev,rep(1:sv, each = q), sep = "_")


  f1 <- Formula(formula)
  mt <- terms(f1, data = data)
  mf1 <- model.frame(f1, data = data)
  y <- as.matrix(model.response(mf1))
  namey <- colnames(y) <- all.vars(f1)[[1]]

  k1 <- NULL
  k2 <- NULL

  if(length(f1)[2L] == 2L){
    x1 <- model.matrix(f1, data = mf1, rhs = 1)
    namesx1 <- colnames(x1)

    x2 <- model.matrix(f1, data = mf1, rhs = 2)
    namesx <- colnames(x2)

    if(any(namesx1 == "(Intercept)") && any(namesx == "(Intercept)"))
      stop("(Intercept) cannot  be specified in both formula")

    k2 <- dim(x2)[2]
    if(any(namesx1 == "(Intercept)")) namesx1[which(namesx1 == "(Intercept)")]  = "Intercept"
    if(any(namesx == "(Intercept)")) namesx[which(namesx == "(Intercept)")] = "Intercept"
    namesxr <- paste(namesx,rep(1:sv,each = k2), sep = "_")

    totc <- sv*k2
    k <- k1 + totc

    xr <- matrix(0, ncol = totc, nrow = n)
    seq_1 <- seq(1, totc, k2)
    seq_2 <- seq(k2, totc,  k2)
    for(i in 1:sv) xr[,seq_1[i]:seq_2[i]] <- x2 * rgm[,i]



    Hmat <- cbind(x1, xr, instr_r)
    Zmat <- cbind(xr, endogvv, x1, endogf)
    colnames(Hmat) <- c(namesx1, namesxr, namesinstr_r)
    colnames(Zmat) <- c(namesxr, namesev_r, namesx1, namesef)


  }
  else{
    x <- model.matrix(f1, data = mf1, rhs = 1)

    k <- dim(x)[2]

    namesx <- colnames(x)
    if(any(namesx == "(Intercept)")) namesx[ which(namesx == "(Intercept)")] = "Intercept"
    namesxr <- paste(namesx,rep(1:sv,each = k), sep = "_")
    totc <- sv*k

    xr <- matrix(0, ncol = totc, nrow = n)
    seq_1 <- seq(1, totc, k)
    seq_2 <- seq(k, totc,  k)

    for(i in 1:sv)     xr[,seq_1[i]:seq_2[i]] <- x * rgm[,i]

    Hmat <- cbind(xr, instr_r)
    Zmat <- cbind(xr, endogvv, endogf)
    colnames(Hmat) <- c(namesxr, namesinstr_r)
    colnames(Zmat) <- c(namesxr, namesev_r, namesef)


  }


  ret <- list(y = y, Hmat = Hmat, Zmat = Zmat)
  return(ret)
}



### S3 methods ----

#' @rdname ivregimes
#' @method coef ivregimes
#' @export
coef.ivregimes <- function(object, ...){
  object[[1]]
}


#' @rdname ivregimes
#' @method vcov ivregimes
#' @import stats
#' @export
vcov.ivregimes <- function(object, ...){
  V <- object[[2]]
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
  # cat("Call:\n")
  # if(is.list((x))) print(x[[1]]$call)
  # else print(x$call)
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
  #   cat("\nCall:\n")
  # if((is.list((x))))  cat(paste(deparse(x[[1]]$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  # else cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

  invisible(x)
}

