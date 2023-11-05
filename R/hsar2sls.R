##### Functions for HSAR model by 2SLS (Based on Chen et al. 2002) ####

#' @title Estimation of HSAR models by 2SLS
#'
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{pdata.frame}.
#' @param listw object. An object of class \code{listw}, \code{matrix}, or \code{Matrix}.
#' @param index index.
#' @param nins numeric. Number of instrument. \code{nins = 2} as default.
#' @param ... additional arguments passed to \code{maxLik}
#' @param x,object an object of class \code{hsar2sls}
#' @param MG logical. If \code{TRUE}, the Mean Group estimator is returned
#' @param digits the number of digits
#' @name hsar2sls
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update))
#' @import stats methods plm spldv
#' @importFrom sphet listw2dgCMatrix
#' @export
hsar2sls <- function(formula,
                     data,
                     listw    = NULL,
                     index    = NULL,
                     nins     = 2,
                     ...)
{
  # ############################ #
  # 1. Initial checks
  # ############################ #

  if (is.null(listw)) stop("listw must be specified")

  # Spatial weight matrix (W): as CsparseMatrix
  if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("Neighbourhood list or listw format unknown")
  if(inherits(listw,"listw"))   W    <- sphet::listw2dgCMatrix(listw)
  if(inherits(listw,"matrix"))  W    <- Matrix(listw)
  if(inherits(listw,"Matrix"))  W    <- listw

  # Check data is pdata.frame
  if (inherits(data, "pdata.frame")  && !is.null(index)) warning("The index argument is ignored because data is a pdata.frame")
  if (!inherits(data, "pdata.frame") && is.null(index)) stop("The data is not pdata.frame and index is needed")
  if (!inherits(data, "pdata.frame"))  data <- pdata.frame(data, index, row.names = FALSE)


  # ############################ #
  # 2. Model Frame
  # ############################ #

  # TODO: add spatially lagged Xs with Formula

  # Model frame
  callT      <- match.call(expand.dots = TRUE)
  callF      <- match.call(expand.dots = FALSE)
  mf         <- callT
  m          <- match(c("formula", "data"), names(mf), 0L)
  mf         <- mf[c(1L, m)]
  mf[[1L]]   <- as.name("model.frame")
  mf$data    <- data
  mf         <- eval(mf, parent.frame())
  nframe     <- length(sys.calls())

  # Get variables and run some checks
  y  <- model.response(mf)
  if (any(is.na(y))) stop("NAs in dependent variable")
  X  <- model.matrix(formula, data = mf, rhs = 1)
  if (any(is.na(X))) stop("NAs in independent variables")

  # IDS
  id   <- attr(data, "index")[[1]]
  tind <- attr(data, "index")[[2]]

  # TODO: Check balanced panel

  # ############################ #
  # 3. 2SLS for each spatial unit
  # ############################ #

  K     <- ncol(X)
  NT    <- nrow(X)
  N     <- length(unique(id))
  TT    <- NT / N

  # TODO: Probably is more efficient to use spreg function

  ## Stack data over i (yt) and then over t (y)
  oo <- order(tind, id)
  y <- y[oo] # All spatial units in each time
  X <- X[oo, ]

  ## Create block diagonal matrix for W
  I_T <- diag(TT)
  Wbd <- kronecker(I_T, W)   # TN * TN

  # Make Wy
  Wy <- Wbd %*% y

  # Make instruments matrix Q
  Xins <- spldv:::make.instruments(Wbd, x = X, q = nins) #TODO: Is this correct?
  H    <- cbind(X, Xins) # NT * p
  # Just linearly independent columns
  H <- H[, qr(H)$pivot[seq_len(qr(H)$rank)]]
  P <- ncol(H)
  if (P < K) stop("Underspecified model")
  if (any(is.na(H))) stop("NAs in the instruments")

  y  <- y[order(oo, id)]
  H  <- H[order(oo, id), ]
  X  <- X[order(oo, id), ]
  Wy <- Wy[order(oo, id), ]
  Z  <- cbind(X, Wy)
  # Generate estimates and standard errors for each spatial unit
  theta.hat <- matrix(NA, nrow = N, ncol = (ncol(X) + 1))
  se        <- matrix(NA, nrow = N, ncol = (ncol(X) + 1))
  se.rob    <- matrix(NA, nrow = N, ncol = (ncol(X) + 1))
  for (i in 1:N){
    anid      <- unique(id)[i]
    theRows   <- which(id == anid)
    Hi        <- H[theRows, ]
    Pi        <- Hi %*% solve(t(Hi)%*% Hi) %*% t(Hi)
    Zi        <- Z[theRows, ]
    yi        <- y[theRows]
    theta.i   <- solve(t(Zi) %*% Pi %*% Zi) %*% t(Zi) %*% Pi %*% yi
    theta.hat[i, ] <- theta.i

    # Standard errors
    r.i      <- yi - crossprod(t(Zi), theta.i)
    sigma2.i <- crossprod(r.i) / (TT - (K + 1))
    tZpZInv  <- solve(t(Zi) %*% Pi %*% Zi)
    V.i      <- drop(sigma2.i) * tZpZInv
    se[i, ]  <- sqrt(diag(V.i))

    # Robust standard errors
    HHinv    <- solve(crossprod(Hi))
    Sigma <- 0
    for (t in 1:TT){
      Sigma       <- Sigma + (r.i[t]^2 * tcrossprod(Hi[t, ]))
    }
    Sigma       <- Sigma / TT
    cheese      <- t(Zi) %*% Hi %*% HHinv %*% Sigma %*% HHinv %*% t(Hi) %*% Zi
    V.ir        <- (TT * (TT / (TT - (K + 1)))) *  tZpZInv %*% cheese %*% tZpZInv
    se.rob[i, ] <- sqrt(diag(V.ir))
  }
  colnames(theta.hat) <- colnames(se) <- colnames(se.rob) <- c(colnames(X), "lambda")
  rownames(theta.hat) <- rownames(se) <- rownames(se.rob) <- unique(id)

  # ############################ #
  # 4. Save results
  # ############################ #
  out <- structure(
    list(
      coefficients = theta.hat,
      y            = y,
      X            = X,
      id           = id,
      tind         = tind,
      N            = N,
      T            = TT,
      se_standard  = se,
      se_robust    = se.rob
    ),
    class = "hsar2sls"
  )
  return(out)
}

## S3 Methods ----

#' @rdname hsar2sls
#' @method summary hsar2sls
#' @export
summary.hsar2sls <- function(object, MG = TRUE, ...){
  if (!MG) stop("Only the standard error of the MG estimator is implemented... We'll work on it")
  N                   <- object$N
  theta               <- object$coefficients
  mean.theta          <- colMeans(theta)
  varMG               <- colSums((theta - repRows(mean.theta, N))^2) / N / (N - 1)
  std.err             <- sqrt(varMG)
  z                   <- mean.theta / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(mean.theta, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  object$MG           <- MG
  class(object)       <- c("summary.hsar2sls", "hsar2sls")
  return(object)
}

#' @rdname hsar2sls
#' @method print summary.hsar2sls
#' @export
print.summary.hsar2sls <- function(x,
                                  digits = max(5, getOption("digits") - 3),
                                  ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                      HSAR by 2SLS \n")
  cat("        ------------------------------------------------------------\n")

  if (x$MG) cat("\nMean Group Estimators :\n")

  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

  cat("Number of obs = ")
  cat(format(x$N * x$T, big.mark = ",", scientific = FALSE))
  cat(sprintf(" (N = %d, T = %d)\n", x$N, x$T))
  invisible(x)
}
