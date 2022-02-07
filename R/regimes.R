


##### Functions for regimes ####
#' Estimation of spatial regime models
#' @name regimes
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{data.frame}.
#' @param splitv variable to split the dataset
#' @param vc   one of ("homoskedastic", "groupwise","heteroskedastic")
#' @param x object for printing
#' @param digits number of digits
#' @param object object for printing
#' @param ... additional arguments to be passed
#'
#' @details
#'
#' The model estimated is:
#'
#' \deqn{
#' y_{ij}= \mathbf{x_{ij,k}}\beta_j + \epsilon
#' }
#' for i=1,..,n representing the sample observations, and j =1,..., J representing
#' the  regimes
#'
#'
#' @author Gianfranco Piras and Mauricio Sarrias
#' @return An object of class ``\code{spregimes}'', a list with elements:
#' \item{coefficients}{the estimated coefficients,}
#' \item{call}{the matched call,}
#' \item{X}{the X matrix}
#' \item{y}{the dependent variable}
#' #' @import Matrix stats spatialreg methods Formula maxLik sphet
#' @export

regimes <- function(formula, data, splitv,
                    vc = c("homoskedastic", "groupwise","heteroskedastic")){

  #Obtain arguments
  vc <- match.arg(vc)
  splitvar <- as.matrix(lm(splitv, data, method="model.frame"))
  sv <- length(unique(splitvar))


  if(vc == "groupwise") res <- groupwise.regimes(formula, data, splitvar)

  else{
    ##do it as a function to accomodate all types of dataset
    nobsg <- table(splitvar)
    csnobsg <- cumsum(nobsg)
    cnt <-  c(1,table(splitvar)+1)[-(sv+1)]
    mt <- terms(formula, data = data)
    mf <- model.frame(formula, data = data)
    data <- data[order(splitvar),]
    y <- model.response(mf)
    x <- model.matrix(mt, mf)
    n <- dim(x)[1]
    k <- dim(x)[2]
    namesx <- colnames(x)
    namesxr <- paste(namesx,rep(1:sv,each = k), sep = "_")
    totr <- sv*k
    xr <- matrix(0, ncol = totr, nrow = n)
    seq_1 <- seq(1, totr, k)
    seq_2 <- seq(k, totr,  k)
    for(i in 1:length(nobsg)) {
      xr[(cnt[i]:csnobsg[i]),seq_1[i]:seq_2[i]] <- x[(cnt[i]):csnobsg[i], ]
    }
    #xr <- as.data.frame(xr)
    names(xr) <- namesxr

  if(vc == "homoskedastic") res <- lm(y~ xr-1, data = data)
  if(vc == "heteroskedastic") stop("hetero not yet implemented")
  }
  class(res) <- "spregimes"
return(res)
  }


groupwise.regimes <- function(formula, data, splitvar){

  dataset <- split.data.frame(data, splitvar)
  res <- lapply(dataset, lm, formula = formula)
  return(res)
}

#' @rdname regimes
#' @method print spregimes
#' @import stats
#' @export
print.spregimes <- function(x,
                         digits = max(3, getOption("digits") - 3),
                         ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print.default(format(drop(coef(x)), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname regimes
#' @method summary spregimes
#' @import stats
#' @export
summary.spregimes <- function(object, ...){
  b                   <- object$coefficients
  std.err             <- sqrt(diag(vcov(object)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.spregime", "spregime")
  return(object)
}

