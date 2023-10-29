##### Functions for HSAR model by QML ####

#' @title Estimation of HSAR models by Quasi-Maximum Likelihood
#'
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{pdata.frame}.
#' @param listw object. An object of class \code{listw}, \code{matrix}, or \code{Matrix}.
#' @param index index.
#' @param gradient logical. Only for testing procedures. Should the analytic gradient be used in the ML optimization procedure? \code{TRUE} as default. If \code{FALSE}, then the numerical gradient is used.
#' @param average  logical. Should the sample log-likelihood function be divided by N?
#' @param init.values if not \code{NULL}, the user must provide a vector of initial parameters for the optimization procedure.
#' @param print.init logical. If \code{TRUE} the initial parameters used in the optimization of the first step are printed.
#' @param otype string. A string indicating whether package \code{maxLik} or \code{optim} is used in for the numerical optimization.
#' @param ... additional arguments passed to \code{maxLik}.
#' @param x,object an object of class \code{hsarML}
#' @param digits the number of digits
#' @name hsarML
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update))
#' @import stats methods maxLik plm
#' @importFrom sphet listw2dgCMatrix
#' @export
hsarML <- function(formula,
                   data,
                   listw       = NULL,
                   index       = NULL,
                   gradient    = TRUE,  # Use the analytical gradient?
                   average     = FALSE, # Use the average log-likelihood function?
                   init.values = NULL,  # vector of starting values
                   print.init  = FALSE,
                   otype       = c("maxLik", "optim"), #type of optimizer
                   ...)
{
  # ############################ #
  # 1. Initial checks
  # ############################ #
  otype        <- match.arg(otype)

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

  # ###################################
  # 2. Model Frame and transform data
  # ###################################

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

  # IDs: spatial unit and time index
  id   <- attr(data, "index")[[1]]
  tind <- attr(data, "index")[[2]]

  # TODO: Check balanced panel

  ## Globals
  K     <- ncol(X)
  NT    <- nrow(X)
  N     <- length(unique(id))
  TT    <- NT / N

  ## Stack data over i (yt) and then over t (y)
  oo <- order(tind, id)
  y <- y[oo]    # All spatial units in each time
  X <- X[oo, ]

  ## Create block diagonal matrix for W and Wy
  I_T <- diag(TT)
  Wbd <- kronecker(I_T, W)   # NT * NT
  Wy  <- crossprod(t(Wbd), y)

  ## Reshape data
  ys  <- matrix(y, nrow = N, ncol = TT)
  Wys <- matrix(Wy, nrow = N, ncol = TT)
  Xs  <- array(0, dim = c(N, TT, K))
  for (k in 1:K){
    Xs[, , k] <- X[ , k]
  }

  # ############################ #
  # 3. Starting Values
  # ############################ #

  if (is.null(init.values)){
    # Lambda
    lambda_init <- cor(drop(Wy), y)
    lambda_init <- rep(lambda_init, N)

    # Beta and sigma
    ols         <- lm(y ~ X - 1)
    beta_init   <- coef(ols)
    beta_init   <- rep(beta_init, N) # coef by individuals
    #sigma2_init <- sum(ols$residuals^2) / ols$df.residual
    sigma2_init <- rep(1, N)

    start <- c(lambda_init, beta_init, sigma2_init) # (N + K*N + N) x 1
    names(start) <- c(paste("lamdba", unique(id), sep = "."),
                      paste(colnames(X), rep(unique(id), each = K), sep = "."),
                      paste("sigma2", unique(id), sep = "."))
  } else {
    start <- init.values
    if (length(start) != N *(K + 2)) stop("Incorrect number of initial parameters")
  }
  if (print.init) print(start)

  sym          <- all(W == t(W))
  omega        <- eigen(W, only.values = TRUE, symmetric = sym)
  lambda_space <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
  if (otype == "maxLik"){
    # Restricted optimization: A %*% theta + B >= 0: Constraint lambda and sigma2
    # maxLik uses a wrapper of optim. We use constrOptim2
    A <- rbind(cbind(diag(1, nrow = N),  matrix(0, nrow = N, ncol = N * K + N)),
               cbind(diag(-1, nrow = N), matrix(0, nrow = N, ncol = N * K + N)),
               cbind(matrix(0, nrow = N, ncol = N),  matrix(0, nrow = N, ncol = N * K), diag(1, nrow = N)))
    B <- c(rep(-1L * (lambda_space[1] + sqrt(.Machine$double.eps)), N),
           rep(lambda_space[2] - sqrt(.Machine$double.eps), N),
           rep(-1L* sqrt(.Machine$double.eps), N))
    # A <- rbind(cbind(diag(1, nrow = N),  matrix(0, nrow = N, ncol = N * K + N)),
    #            cbind(diag(-1, nrow = N), matrix(0, nrow = N, ncol = N * K + N))
    #            )
    # B <- c(rep(-1* (lambda_space[1] + .Machine$double.eps), N),
    #        rep(lambda_space[2]      - .Machine$double.eps,  N))
    callT$constraints <- list(ineqA = A, ineqB = B)
  } else {
    LB <- c(rep(lambda_space[1] + sqrt(.Machine$double.eps), N),
            rep(-Inf, N * K),
            rep(.Machine$double.eps, N))
    UB <- c(rep(lambda_space[2] - sqrt(.Machine$double.eps), N),
            rep(Inf, N * K),
            rep(Inf, N))
    # LB <- c(rep(-0.995, N),
    #        rep(-Inf, N * K),
    #        rep(0.010, N))
    # UB <- c(rep(0.995, N),
    #        rep(Inf, N * K),
    #        rep(Inf, N))
  }

  # ################
  # 4. Optimization
  # ################

  if (otype == "maxLik"){
    # Optimization default controls if not added by user
    if (is.null(callT$method)) callT$method  <- 'bfgs'
    #if (is.null(callT$rho)) callT$rho  <- 'bfgs'
    if (is.null(callT$iterlim)) callT$iterlim <- 100000
    callT$finalHessian <- FALSE # We do not require the Hessian. This speeds the optimization procedure.

    opt <- callT
    m <- match(c('method', 'print.level', 'iterlim',
                 'tol', 'ftol', 'steptol', 'fixed', 'constraints',
                 'control', 'finalHessian', 'reltol', 'rho', 'outer.iterations', 'outer.eps'),
               names(opt), 0L)
    opt <- opt[c(1L, m)]
    opt$start     <- start
    opt[[1]]      <- as.name('maxLik')
    opt$logLik    <- as.name('mlfunc')
    opt$average   <- as.name('average')
    opt$gradient  <- gradient
    opt$post      <- FALSE
    opt[c('y', 'W', 'X', 'id', 'Wy', 'oo')] <- list(as.name('ys'),
                                                    as.name('W'),
                                                    as.name('Xs'),
                                                    as.name('id'),
                                                    as.name('Wys'),
                                                    as.name('oo'))
    out <- eval(opt, sys.frame(which = nframe))
  } else  {
    # Optim uses "L-BFGS-B"
    if (is.null(callT$method)) callT$method  <- 'L-BFGS-B'
    if (is.null(callT$maxit))  callT$maxit   <- 100000
    opt <- callT
    m1 <- match(c('method', 'hessian'), names(opt), 0L)
    m2 <- match(c('reltol', 'abstol', 'REPORT', 'maxit', 'trace', 'pgtol'),
                names(opt), 0L)
    control <- as.list(opt[c(m2)])
    opt     <- opt[c(1L, m1)]
    opt$par      <- start
    opt[[1]]     <- as.name('optim')
    opt$fn       <- as.name('ml_optim')
    opt$gr       <- as.name('gr_optim')
    opt$lower    <- LB
    opt$upper    <- UB
    opt$average  <- average
    opt$control  <- control
    opt[c('y', 'W', 'X', 'id', 'Wy', 'oo')] <- list(as.name('ys'),
                                                    as.name('W'),
                                                    as.name('Xs'),
                                                    as.name('id'),
                                                    as.name('Wys'),
                                                    as.name('oo'))
    out <- eval(opt, sys.frame(which = nframe))
  }


  # ############################ #
  # 5. Save results
  # ############################ #

  # Obtain Hessian and Git
  theta.hat       <- if (otype == "maxLik") coef(out) else out$par
  opt[[1]]        <- as.name('mlfunc')
  names(opt)[[2]] <- 'theta'
  opt$post        <- TRUE
  opt$gradient    <- FALSE
  opt[[2]]        <- theta.hat
  again           <- eval(opt, sys.frame(which = nframe))
  hessian         <- attr(again, 'hessian')
  Git             <- attr(again, 'Git')

  # # Classical standard errors
  InvH            <- solve(hessian)
  V               <- InvH / TT # Note that V_hat = V_asym / TT
  se              <- sqrt(diag(V))

  # Robust standard errors
  J      <- (Git %*% t(Git)) / TT
  V.rob  <- (InvH %*% J %*% InvH) / TT
  se.rob <- sqrt(diag(V.rob))

  #Save standard errors
  se.lambda  <- se[1:N] # N x 1
  se.betas   <- se[(N + 1):(N + (K * N))] # KN x 1
  se.sigmas2 <- se[(N + (K * N) + 1):(N + (N * K) + N)] # N x 1
  se.Betas   <- matrix(se.betas, K, N)
  se.theta   <- cbind(se.lambda, t(se.Betas), se.sigmas2)
  colnames(se.theta)  <- c("lambda", colnames(X), "sigma2")
  rownames(se.theta)  <- unique(id)
  #se.theta <- NULL

  #Save robust standard errors
  ser.lambda  <- se.rob[1:N] # N x 1
  ser.betas   <- se.rob[(N + 1):(N + (K * N))] # KN x 1
  ser.sigmas2 <- se.rob[(N + (K * N) + 1):(N + (N * K) + N)] # N x 1
  ser.Betas   <- matrix(ser.betas, K, N)
  ser.theta   <- cbind(ser.lambda, t(ser.Betas), ser.sigmas2)
  colnames(ser.theta)  <- c("lambda", colnames(X), "sigma2")
  rownames(ser.theta)  <- unique(id)
  #ser.theta <- NULL

  # Theta: (lambdas, betas, sigmas)
  lambda  <- theta.hat[1:N] # N x 1
  betas   <- theta.hat[(N + 1):(N + (K * N))] # KN x 1
  sigmas2 <- theta.hat[(N + (K * N) + 1):(N + (N * K) + N)] # N x 1
  Betas   <- matrix(betas, K, N)
  theta   <- cbind(lambda, t(Betas), sigmas2)
  colnames(theta)  <- c("lambda", colnames(X), "sigma2")
  rownames(theta)  <- unique(id)

  maxlik <- out
  out <- structure(
    list(
      coefficients = theta,
      y            = y,
      X            = X,
      id           = id,
      tind         = tind,
      N            = N,
      T            = TT,
      hessian      = hessian,
      Git          = Git,
      se_standard  = se.theta,
      se_robust    = ser.theta,
      maxlik       = maxlik
    ),
    class = "hsarML"
  )
  return(out)
}


## S3 Methods ----
#' @rdname hsarML
#' @export
coef.hsarML <- function(object, ...){
  object$coefficients
}

#' @rdname hsarML
#' @method summary hsarML
#' @export
summary.hsarML <- function(object, MG = TRUE, ...){
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
  class(object)       <- c("summary.hsarML", "hsarML")
  return(object)
}

#' @rdname hsarML
#' @method print summary.hsarML
#' @export
print.summary.hsarML <- function(x,
                                 digits = max(5, getOption("digits") - 3),
                                 ...)
{
  cat("        ------------------------------------------------------------\n")
  cat("                      HSAR by ML \n")
  cat("        ------------------------------------------------------------\n")

  if (x$MG) cat("\nMean Group Estimators :\n")

  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

  cat("Number of obs = ")
  cat(format(x$N * x$T, big.mark = ",", scientific = FALSE))
  cat(sprintf(" (N = %d, T = %d)\n", x$N, x$T))
  invisible(x)
}

## Additional functions ----
repCols <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

repRows <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}


## QML functions ----
mlfunc <- function(theta, y, Wy, W, X, id, oo,
                   gradient  = TRUE,
                   post      = FALSE,
                   average   = FALSE,
                   ...){
  # Globals
  K     <- dim(X)[3]
  TT    <- dim(X)[2]
  N     <- length(unique(id))
  NT    <- TT * N

  # Theta: the order of the parameters is : lambdas, betas, sigmas
  lambdas <- theta[1:N]                                 # N x 1
  betas   <- theta[(N + 1):(N + (K * N))]               # KN x 1
  sigmas2 <- theta[(N + (K * N) + 1):(N + (N * K) + N)] # N x 1
  Betas   <- matrix(betas, K, N)                        # K x N

  ## Create XB
  Xb  <- matrix(0, N, TT)
  for (k in 1:K){
    Xb_k <- X[, , k] * rep(Betas[k, ], TT)  #N x T
    Xb  <- Xb + Xb_k
  }
  res.it <- y - Wy * rep(lambdas, TT) - Xb
  res2.i <- rowSums(res.it^2)

  # Determinant
  # TODO: Add Ord's approximation
  A    <- diag(N) - crossprod(t(diag(lambdas)), W)
  detA <- det(A)
  if(detA <= 0){
    stop("`I -  Lambda W` is not invertible")
  }
  ## See Equation 46 in Aquaro et al (2015)
  p1     <- - log(2 * pi) * NT / 2
  p2     <- - (TT / 2) * sum(log(sigmas2))
  p3     <- TT * log(detA)
  p4     <- - (1 / 2) * sum(res2.i / sigmas2)
  LL     <- if (average) (p1 + p2 + p3 + p4) / TT else (p1 + p2 + p3 + p4)


  if (gradient){
    # dl/dlambda
    G       <- W %*% solve(A)
    resWy.i <- rowSums(res.it * Wy)
    dlambda <- - TT * diag(t(G)) + (resWy.i / sigmas2)
    # dl/db
    dbetas <- matrix(NA, N, K)
    for (k in 1:K){
      dbetak <- rowSums(X[, , k] * res.it) / sigmas2
      dbetas[, k] <- dbetak
    }
    # dl/s2
    dsig2   <- - (TT / (2 * sigmas2)) +  res2.i / (2 * sigmas2^2)
    g.bar   <- c(dlambda, as.numeric(t(dbetas)), dsig2)
    attr(LL, 'gradient') <- if (average) g.bar / TT else g.bar
  }
  if (post){
    # Create Git which is N(K + 2) x T
    # Create Hessian: -(1/T) H
    # dl/dlambda
    G       <- W %*% solve(A)
    wy2.i   <- rowSums(Wy^2)
    #wy2.i   <- apply(as.matrix(Wy.s^2), 2, tapply, id, sum)
    resWy.i  <-  rowSums(res.it * Wy)
    #resWy.i <- apply(as.matrix(res.it * Wy[order(oo, id)]), 2, tapply, id, sum)   # N x 1

    H11   <- G * t(G) + diag(drop(wy2.i/ sigmas2), nrow = N) / TT                                      # N x N
    H13   <- diag(drop(resWy.i / sigmas2^2), nrow = N) / TT                                             # N x N
    H33   <- diag((- 1 / (2 * sigmas2^2)), nrow = N) + diag(drop(res2.i / (sigmas2^3)), nrow = N) / TT  # N x N
    H12   <- matrix(0, N , N * K)
    H22   <- matrix(0, N * K, N * K)
    H23   <- matrix(0, N * K, N)

    dlambda.it <- matrix(NA, N, TT)
    dbetas.it  <- matrix(NA, N * K, TT)
    dsigmas.it <- matrix(NA, N , TT)
    for (i in 1:N){
      #anid             <- unique(id)[i]
      #theRows          <- which(id == anid)
      cind             <- c((i - 1) * K + 1):(i * K)
      Wy.i             <- Wy[i, ]
      X.i              <- X[i, , ]   # T x K
      s.i              <- sigmas2[i]
      r.i              <- res.it[i, ]
      H12[i, cind]     <- t(Wy.i) %*% (X.i / s.i) / TT  # 1 x K
      H22[cind, cind]  <- t(X.i) %*% X.i / s.i  /  TT   # 1 x K
      H23[cind, i]     <- t(X.i) %*%  r.i / s.i^2 / TT  # K x 1
      dlambda.it[i, ]  <- - diag(G)[i] + (r.i * Wy.i / s.i)
      dbetas.it[cind, ] <- t(r.i * X.i / s.i)
      dsigmas.it[i, ]   <- - (1 / (2 * s.i)) +  r.i^2 / (2 * s.i^2)
    }
    row_1   <- cbind(H11, H12, H13)
    row_2   <- cbind(t(H12), H22, H23)
    row_3   <- cbind(t(H13), t(H23), H33)
    hessian <- rbind(row_1, row_2, row_3)
    attr(LL, 'hessian') <- hessian
    Git   <- rbind(dlambda.it, dbetas.it, dsigmas.it)
    attr(LL, 'Git') <- Git
  }
  return(LL)
}

ml_optim <- function(theta, y, W, Wy, X, id, oo,
                     average   = FALSE, ...){
  K     <- dim(X)[3]
  TT    <- dim(X)[2]
  N     <- length(unique(id))
  NT    <- TT * N

  # Theta: the order of the parameters is : lambdas, betas, sigmas
  lambdas <- theta[1:N]                                 # N x 1
  betas   <- theta[(N + 1):(N + (K * N))]               # KN x 1
  sigmas2 <- theta[(N + (K * N) + 1):(N + (N * K) + N)] # N x 1
  Betas   <- matrix(betas, K, N)                        # K x N

  ## Create XB
  Xb  <- matrix(0, N, TT)
  for (k in 1:K){
    Xb_k <- X[, , k] * rep(Betas[k, ], TT)  #N x T
    Xb  <- Xb + Xb_k
  }
  res.it <- y - Wy * rep(lambdas, TT) - Xb
  res2.i <- rowSums(res.it^2)

  # Determinant
  # TODO: Add Ord's approximation
  A    <- diag(N) - crossprod(t(diag(lambdas)), W)
  detA <- det(A)
  if(detA <= 0){
    stop("`I -  Lambda W` is not invertible")
  }
  ## See Equation 46 in Aquaro et al (2015)
  p1     <- - log(2 * pi) * NT / 2
  p2     <- - (TT / 2) * sum(log(sigmas2))
  p3     <- TT * log(detA)
  p4     <- - (1 / 2) * sum(res2.i / sigmas2)
  LL     <- if (average) (p1 + p2 + p3 + p4) / TT else (p1 + p2 + p3 + p4)

  return(-1 * LL)
}

gr_optim <- function(theta, y, W, Wy, X, id, oo,
                     average   = FALSE,
                     ...){
  # Globals
  K     <- dim(X)[3]
  TT    <- dim(X)[2]
  N     <- length(unique(id))
  NT    <- TT * N

  # Theta: the order of the parameters is : lambdas, betas, sigmas
  lambdas <- theta[1:N]                                 # N x 1
  betas   <- theta[(N + 1):(N + (K * N))]               # KN x 1
  sigmas2 <- theta[(N + (K * N) + 1):(N + (N * K) + N)] # N x 1
  Betas   <- matrix(betas, K, N)                        # K x N

  ## Create XB
  Xb  <- matrix(0, N, TT)
  for (k in 1:K){
    Xb_k <- X[, , k] * rep(Betas[k, ], TT)  #N x T
    Xb  <- Xb + Xb_k
  }
  res.it <- y - Wy * rep(lambdas, TT) - Xb
  res2.i <- rowSums(res.it^2)

  # Determinant
  # TODO: Add Ord's approximation
  A    <- diag(N) - crossprod(t(diag(lambdas)), W)
  detA <- det(A)
  if(detA <= 0){
    stop("`I -  Lambda W` is not invertible")
  }
  ## See Equation 46 in Aquaro et al (2015)
  p1     <- - log(2 * pi) * NT / 2
  p2     <- - (TT / 2) * sum(log(sigmas2))
  p3     <- TT * log(detA)
  p4     <- - (1 / 2) * sum(res2.i / sigmas2)
  LL     <- if (average) (p1 + p2 + p3 + p4) / TT else (p1 + p2 + p3 + p4)


  # dl/dlambda
  G       <- W %*% solve(A)
  resWy.i <- rowSums(res.it * Wy)
  dlambda <- - TT * diag(t(G)) + (resWy.i / sigmas2)
  # dl/db
  dbetas <- matrix(NA, N, K)
  for (k in 1:K){
    dbetak <- rowSums(X[, , k] * res.it) / sigmas2
    dbetas[, k] <- dbetak
  }
  # dl/s2
  dsig2   <- - (TT / (2 * sigmas2)) +  res2.i / (2 * sigmas2^2)
  g.bar   <- c(dlambda, as.numeric(t(dbetas)), dsig2)
  g.bar   <- if (average) g.bar / TT else g.bar

  return(g.bar * -1)
}


