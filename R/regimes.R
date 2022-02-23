


##### Functions for regimes ####
#' Estimation of spatial regime models
#' @name regimes
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{data.frame}.
#' @param rgv variable to identify the regimes
#' @param vc   one of ("homoskedastic", "groupwise")
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
#' @return An object of class ``\code{lm}'', a list with elements:
#' @import Formula sphet stats
#' @export

regimes <- function(formula, data, rgv = NULL,
                    vc = c("homoskedastic", "groupwise")){



  if(is.null(rgv)) stop("regimes variable not specified")
  if(class(rgv) != "formula") stop("regimes variable has to be a formula")

  #Obtain arguments
  vc                <- match.arg(vc)

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

  #class(res) <- "spregimes"
return(res)
  }


groupwise.regimes <- function(formula, data, k, k1, k2, rgm, sv){

  if (is.null(k2)) df <- k
  else df <- k1 + 2*k2


  fs <- lm(formula, data)
  uhat <- residuals(fs)
  nobsg <- colSums(rgm)

  omega <- vector("numeric", length = nrow(data))
  for (i in 1: sv) omega[which(as.numeric(rgm[,i])==1)] <-  1/(crossprod(residuals(fs)[which(as.numeric(rgm[,i])==1)])/(nobsg[i]-df))

  data$omega <- omega
  res <-  lm(formula, data, weights = omega)

  return(res)
}

homoskedastic.regimes <- function(formula, dataset){
  res <- lm(formula, dataset)
  return(res)
}
data.prep.regimes <- function(formula, data, rgv){

  splitvar         <- as.matrix(lm(rgv, data, method="model.frame"))
  sv               <- length(unique(splitvar))
  svm              <- as.numeric(unique(splitvar))
  rgm              <-  matrix(,nrow = nrow(data), ncol = 0)
  for(i in svm)    rgm <- cbind(rgm, ifelse(splitvar ==  i, 1, 0))



  f1 <- Formula(formula)
  mt <- terms(f1, data = data)
  mf1 <- model.frame(f1, data = data)
  y <- model.response(mf1)

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
    namesxr <- paste(namesx,rep(1:sv,each = k2), sep = "_")

    totc <- sv*k2
    k <- k1 + totc

    xr <- matrix(0, ncol = totc, nrow = n)
    seq_1 <- seq(1, totc, k2)
    seq_2 <- seq(k2, totc,  k2)
    for(i in 1:sv) xr[,seq_1[i]:seq_2[i]] <- x2 * rgm[,i]

    data <- data.frame(y, x1,xr)
    colnames(data) <- c("y", namesx1, namesxr)
    form <- as.formula(paste("y ~", paste(names(data)[2:(k+1)], collapse =  " + "), "-1"))

  }
  else{
    x <- model.matrix(f1, data = mf1, rhs = 1)
    n <- dim(x)[1]
    k <- dim(x)[2]
    namesx <- colnames(x)
    if(any(namesx == "(Intercept)")) namesx[ which(namesx == "(Intercept)")] = "Intercept"
    namesxr <- paste(namesx,rep(1:sv,each = k), sep = "_")
    totc <- sv*k

    xr <- matrix(0, ncol = totc, nrow = n)
    seq_1 <- seq(1, totc, k)
    seq_2 <- seq(k, totc,  k)

    for(i in 1:sv)     xr[,seq_1[i]:seq_2[i]] <- x * rgm[,i]
    data <- data.frame(y,xr)

    colnames(data) <- c("y" ,namesxr)

    if("Intercept_1" %in% namesxr)  form <- as.formula(paste("y ~", paste(names(data)[2:(totc+1)], collapse =  " + "), "-1"))
    else form <- as.formula(paste("y ~", paste(names(data)[2:(totc+1)], collapse =  " + ")))

  }


  ret <- list(data, form, k, k1, k2, rgm, sv)
return(ret)
}
