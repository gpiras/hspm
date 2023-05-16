##### Functions for regimes ####
#' @title Estimation of regime models with endogenous variables
#' @description The function \code{ivregimes} deals with
#' the estimation of regime models.
#' Most of the times the variable identifying the regimes
#' reveals some spatial aspects of the data (e.g., administrative boundaries).
#' The model includes exogenous as well as endogenous
#' variables among the regressors.

#' @name ivregimes
#' @param formula a symbolic description of the model of the form \code{y ~ x_f | x_v | h_f | h_v} where \code{y} is the dependent variable, \code{x_f} are the regressors that do not vary by regimes,  \code{x_v} are the regressors that vary by regimes, \code{h_f} are the fixed instruments and \code{h_v} are the instruments that vary by regimes.
#' @param data the data of class \code{data.frame}.
#' @param rgv an object of class \code{formula} to identify the regime variables
#' @param vc   one of \code{c("homoskedastic", "robust", "OGMM")}. If \code{"OGMM"} an optimal weighted GMM is used to estimate the VC matrix.
#'
#'
#' @details
#'
#' The basic (non spatial) model with endogenous variables
#' can be written in a general way as:
#'  \deqn{
#' y
#' =
#' \begin{bmatrix}
#' X_1& 0 \\
#' 0 & X_2 \\
#' \end{bmatrix}
#' \begin{bmatrix}
#' \beta_1 \\
#' \beta_2 \\
#' \end{bmatrix}
#' + X\beta +
#' \begin{bmatrix}
#' Y_1& 0 \\
#' 0 & Y_2 \\
#' \end{bmatrix}
#' \begin{bmatrix}
#' \pi_1 \\
#' \pi_2 \\
#' \end{bmatrix}
#' + Y\pi +
#'  \varepsilon
#' }
#' where  \eqn{y = [y_1^\prime,y_2^\prime]^\prime},
#' and the \eqn{n_1 \times 1} vector \eqn{y_1} contains the observations
#' on the dependent variable for the first regime,
#' and the \eqn{n_2 \times 1} vector \eqn{y_2} (with \eqn{n_1 + n_2 = n})
#' contains the observations on the dependent variable for the second regime.
#' The \eqn{n_1 \times k} matrix \eqn{X_1} and the \eqn{n_2 \times k}
#' matrix \eqn{X_2} are blocks of a block diagonal matrix,
#' the vectors of parameters  \eqn{\beta_1} and \eqn{\beta_2} have
#' dimension \eqn{k_1 \times 1} and \eqn{k_2 \times 1}, respectively,
#' \eqn{X} is the \eqn{n \times p} matrix of regressors that do not vary by regime,
#' \eqn{\beta}  is a \eqn{p\times 1} vector of parameters.
#' The three matrices \eqn{Y_1} (\eqn{n_1 \times q}),
#' \eqn{Y_2} (\eqn{n_2 \times q}) and \eqn{Y} (\eqn{n \times r})
#' with corresponding vectors of parameters \eqn{\pi_1}, \eqn{\pi_2} and \eqn{\pi},
#' contain the endogenous variables.
#' Finally, \eqn{\varepsilon = [\varepsilon_1^\prime,\varepsilon_2^\prime]^\prime}
#' is the \eqn{n\times 1} vector of innovations.
#' The model is estimated by two stage least square.
#' In particular:
#' \itemize{
#' \item If \code{vc = "homoskedastic"},
#' the variance-covariance matrix is estimated by \eqn{\sigma^2(\hat Z^\prime \hat Z)^{-1}},
#' where \eqn{\hat Z= PZ},  \eqn{P= H(H^\prime H)^{-1}H^\prime}, \eqn{H} is the matrix of instruments,
#' and \eqn{Z} is the matrix of all exogenous and endogenous variables in the model.
#'
#' \item If \code{vc = "robust"}, the variance-covariance matrix is estimated by
#' \eqn{(\hat Z^\prime \hat Z)^{-1}(\hat Z^\prime \hat\Sigma \hat Z) (\hat Z^\prime \hat Z)^{-1}},
#' where \eqn{\hat\Sigma} is a diagonal matrix with diagonal elements \eqn{\hat\sigma_i},
#' for \eqn{i=1,...,n}.
#' \item Finally, if \code{vc = "OGMM"}, the model is estimated in two steps.
#' In the first step, the model is estimated by 2SLS yielding
#' the residuals \eqn{\hat \varepsilon}.
#' With the residuals, the diagonal matrix \eqn{\hat \Sigma} is estimated and is
#' used to construct the matrix \eqn{\hat S = H^\prime \hat \Sigma H}.
#' Then \eqn{\eta_{OWGMM}=(Z^\prime H\hat S^{-1}H^\prime Z)^{-1}Z^\prime H\hat S^{-1}H^\prime y}, where \eqn{\eta_{OWGMM}}
#' is the vector of all the parameters in the model,
#' The variance-covariance matrix is: \eqn{n(Z^\prime H\hat S^{-1}H^\prime Z)^{-1}}.}

#' @examples
#' data("natreg")
#' form   <- HR90  ~ 0 | MA90 + PS90 + RD90 + UE90 | 0 | MA90 + PS90 + RD90 + FH90 + FP89 + GI89
#' split  <- ~ REGIONS
#' mod <- ivregimes(formula = form, data = natreg, rgv = split, vc = "robust")
#' summary(mod)
#' mod1 <- ivregimes(formula = form, data = natreg, rgv = split, vc = "OGMM")
#' summary(mod1)
#' form1   <- HR90  ~ MA90 + PS90 |  RD90 + UE90 -1 | MA90 + PS90 | RD90 + FH90 + FP89 + GI89 -1
#' mod2 <- ivregimes(formula = form1, data = natreg, rgv = split, vc = "homoskedastic")
#' summary(mod2)
#'
#' @author Gianfranco Piras and Mauricio Sarrias
#' @return An object of class \code{ivregimes}. A \code{list} of five elements. The first element of the list contains the estimation results. The other elements are needed for printing the results.
#' @import Formula sphet stats
#' @export

ivregimes <- function(formula, data, rgv = NULL,
                      vc = c("homoskedastic", "robust", "OGMM")){

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
  ct       <- intro[[7]][[6]]

  res <- tsls_regimes(y, Hmat, Zmat, vc)

  res <- Matchgroups(res, ct)

  colnames.end.v <- Matchnames(colnames.end.v, ct)
  colnames.instr <- Matchnames(colnames.instr, ct)


  res <- list(res, cl, colnames.end.f, colnames.end.v, colnames.instr)
  class(res) <- c("spregimes","ivregimes")
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
  yp <- array(yp, dim = c(length(yp),1),
              dimnames = list(seq(1, length(yp)), ""))

  e  <- as.matrix(y) - yp
  e <- array(e, dim = c(length(e),1),
             dimnames = list(seq(1, length(e)), ""))


  ##simple
  if(vc == "homoskedastic"){

    cpe <- crossprod(e)
    vcmatrix <- (as.numeric(cpe) /df) * ZpZpi
    results <- list(coefficients = betaiv,
                    var = vcmatrix,
                    residuals = e,  X = Zmat,
                    y = y, yp = yp)

  }


  if(vc == "robust"){
    e2 <- e^2
    ZoZ <- crossprod(Zp, (as.matrix(Zp) * as.numeric(e2)))
    vcmatrix <-  ZpZpi %*% ZoZ %*% ZpZpi * (n/df)
    results <- list(coefficients = betaiv,
                    var = vcmatrix,
                    residuals = e,  X = Zmat,
                    y = y, yp = yp)

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
    results <- list(coefficients = b_OWGMM,
                    var = vcmatrix,
                    residuals = e,  X = Zmat,
                    y = y, yp = yp)

  }
  return(results)
}


iv.data.prep.regimes <- function(formula, data, rgv){
  n                <- dim(data)[1]
  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- sort(as.numeric(unique(splitvar)), decreasing = F)
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))


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
  namesxv <- paste(namesxv, rep(1:sv, each = k2), sep = "_")
  colnames(XV) <- namesxv
  namesZ <- c(namesxf, namesxv)


  ## extract instrument
  Zf <- model.matrix(F1, data = mf, rhs = 3, drop = FALSE)
  nameszf <- colnames(Zf)
  Zv <- model.matrix(F1, data = mf, rhs = 4, drop = FALSE)
  nameszv <- colnames(Zv)

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
  nameszV <- paste(nameszv,rep(1:sv,each = k3), sep = "_")
  namesH <- c(nameszf, nameszV)




  ## if instruments check identification
  end.f <- Xf[, !(colnames(Xf) %in% colnames(Zf)), drop = FALSE]
  end.v <- Xv[, !(colnames(Xv) %in% colnames(Zv)), drop = FALSE]

  colnames.end.f <- colnames(end.f)
  colnames.end.v <- colnames(end.v)
  colnames.end.v <- paste(rep(colnames.end.v,  each = sv), "_", 1:sv, sep = "")

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

  ret <- list(y = y, Hmat = Hmat, Zmat = Zmat, colnames.end.f, colnames.end.v, colinst, l.split)
  return(ret)
}




