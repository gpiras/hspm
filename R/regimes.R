


##### Functions for regimes ####
#' @title Estimation of regimes models
#' @description The function \code{regimes} deals with
#' the estimation of regime models.
#' Most of the times the variable identifying the regimes
#' reveals some spatial aspects of the data (e.g., administrative boundaries).
#'
#'
#' @name regimes
#' @param formula a symbolic description of the model of the form \code{y ~ x_f | x_v} where \code{y} is the dependent variable, \code{x_f} are the regressors that do not vary by regimes and  \code{x_v} are the regressors that vary by regimes
#' @param data the data of class \code{data.frame}.
#' @param rgv an object of class \code{formula} to identify the regime variables
#' @param vc one of \code{c("homoskedastic", "groupwise")}. If \code{groupwise}, the model VC matrix is estimated by weighted least square.
#' @details
#' For convenience and without loss of generality,
#' we assume the presence of only two regimes.
#' In this case,
#' the basic (non-spatial) is:

#'\deqn{
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
#' \varepsilon
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
#' \eqn{\beta}  is a \eqn{p\times 1} vector of parameters
#' and \eqn{\varepsilon = [\varepsilon_1^\prime,\varepsilon_2^\prime]^\prime}
#' is the \eqn{n\times 1} vector of innovations.
#' \itemize{
#' \item If \code{vc = "homoskedastic"}, the model is estimated by OLS.
#' \item If \code{vc = "groupwise"}, the model is estimated in two steps.
#' In the first step, the model is estimated by OLS. In the second step, the
#' inverse of the (groupwise) residuals from the first step are employed
#' as weights in a  weighted least square procedure.}
#'
#' @examples
#' data("baltim")
#' form   <- PRICE  ~ NROOM + NBATH + PATIO + FIREPL + AC + GAR + AGE + LOTSZ + SQFT
#' split  <- ~ CITCOU
#' mod <- regimes(formula = form, data = baltim, rgv = split, vc = "groupwise")
#' summary(mod)
#' form <- PRICE  ~ AC + AGE + NROOM + PATIO + FIREPL + SQFT | NBATH + GAR + LOTSZ - 1
#' mod <- regimes(form, baltim, split, vc = "homoskedastic")
#' summary(mod)
#'
#'
#' @author Gianfranco Piras and Mauricio Sarrias
#' @return An object of class \code{lm} and \code{spregimes}.
#' @import Formula sphet stats spdep
#' @export

regimes <- function(formula, data, rgv = NULL,
                    vc = c("homoskedastic", "groupwise")){



  if(is.null(rgv)) stop("regimes variable not specified")
  if(!inherits(rgv, "formula")) stop("regimes variable has to be a formula")

  #Obtain arguments
  vc                <- match.arg(vc)
  cl                <- match.call()

  #process the data
  intro    <- data.prep.regimes(formula = formula, data = data, rgv = rgv)

  #split the object list intro:
  # data
  # formula
  # k tot number of variables in model
  # totr tot of variables varying by regime

  dataset  <- intro[[1]]
  form     <- intro[[2]]
  k        <- intro[[3]]
  k1       <- intro[[4]]
  k2       <- intro[[5]]
  rgm      <- intro[[6]]
  sv       <- intro[[7]]

if (vc == "groupwise")  res <- groupwise.regimes(form, dataset, k, k1, k2, rgm, sv)

if (vc == "homoskedastic")  res <- lm(form, dataset)

 res <- list(res, cl)

class(res) <- c("spregimes", "lm")
return(res)
  }


groupwise.regimes <- function(formula, data, k, k1, k2, rgm, sv){

  if (is.null(k2)) df <- k
  else df <- k1 + 2*k2


  fs <- lm(formula, data)
  nobsg <- colSums(rgm)

  omega <- vector("numeric", length = nrow(data))
  for (i in 1: sv) omega[which(as.numeric(rgm[,i])==1)] <-  1/(crossprod(residuals(fs)[which(as.numeric(rgm[,i])==1)])/(nobsg[i]-df))

  data$omega <- omega

  res <-  lm(formula, data, weights = omega)
  # uhat <- vector("numeric", length = nrow(data))
  # for (i in 1: sv) uhat[which(as.numeric(rgm[,i])==1)] <-  residuals(lm(formula, data[which(as.numeric(rgm[,i])==1),] ))
  # fs$residuals <- uhat

  return(res)
}

###homoskedastic.regimes <- function(formula, dataset){

#  res <- lm(formula, dataset)
#  return(res)

  #}
##
data.prep.regimes <- function(formula, data, rgv){

  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- sort(as.numeric(unique(splitvar)), decreasing = F)
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))



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
    if(any(namesx1 == "(Intercept)") && any(colnames(x2) == "(Intercept)"))
                stop("(Intercept) cannot  be specified in both formula")
    n <- dim(x1)[1]
    k1 <- dim(x1)[2]
    k2 <- dim(x2)[2]
   if(any(namesx1 == "(Intercept)")) namesx1[which(namesx1 == "(Intercept)")]  = "Intercept"
   if(any(namesx == "(Intercept)")) namesx[which(namesx == "(Intercept)")] = "Intercept"

    namesxr <- paste(namesx, rep(svm, each = k2), sep = "_")
    totc <- sv*k2
    k <- k1 + totc

    xr <- matrix(0, ncol = totc, nrow = n)
    seq_1 <- seq(1, totc, k2)
    seq_2 <- seq(k2, totc,  k2)
    for(i in 1:sv) xr[,seq_1[i]:seq_2[i]] <- x2 * rgm[,i]

    data <- data.frame(y, x1,xr)
    colnames(data) <- c(namey, namesx1, namesxr)
    form <- as.formula(paste(namey," ~", paste(names(data)[2:(k+1)], collapse =  " + "), "-1"))

  }
  else{
    x <- model.matrix(f1, data = mf1, rhs = 1)
    n <- dim(x)[1]
    k <- dim(x)[2]
    namesx <- colnames(x)

    Int    <- FALSE
    if(any(namesx == "(Intercept)")) {
      Int <- TRUE
      namesx[ which(namesx == "(Intercept)")] = "Intercept"
      Int <- TRUE
      }

    namesxr <- paste(namesx, rep(svm, each = k), sep = "_")

    totc <- sv*k

    xr <- matrix(0, ncol = totc, nrow = n)
    seq_1 <- seq(1, totc, k)
    seq_2 <- seq(k, totc,  k)

    for(i in 1:sv)     xr[,seq_1[i]:seq_2[i]] <- x * rgm[,i]
    data <- data.frame(y,xr)

    colnames(data) <- c(namey ,namesxr)

    if(isTRUE(Int))  form <- as.formula(paste(namey,"~", paste(names(data)[2:(totc+1)], collapse =  " + "), "-1"))
    else form <- as.formula(paste(namey,"~", paste(names(data)[2:(totc+1)], collapse =  " + ")))

  }


  ret <- list(data, form, k, k1, k2, rgm, sv)
return(ret)
}



