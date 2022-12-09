
### S3 methods ----
#' @rdname spregimes
#' @method coef spregimes
#' @export
coef.spregimes <- function(object, ...){
  coeff<- object[[1]][[1]]
  if("lm" %in% class(object)) coeff <- coefficients(object[[1]])
  return(coeff)
}


#' @rdname spregimes
#' @method vcov spregimes
#' @import stats
#' @export
vcov.spregimes <- function(object, ...){
  V <- object[[1]][[2]]
  if("lm" %in% class(object))   V <- vcov(object[[1]])
  return(V)
}




#' @rdname spregimes
#' @method print spregimes
#' @import stats
#' @export
print.spregimes <- function(x,
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



#' @rdname spregimes
#' @method summary spregimes
#' @import stats
#' @export
summary.spregimes <- function(object, ...){
  b                   <- coef(object)
  std.err             <- sqrt(diag(vcov(object)))
  z                   <- b / std.err
  p                   <- 2 * (1 - pnorm(abs(z)))
  CoefTable           <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object)       <- c("summary.spregimes", class(object)[2])
  return(object)
}


#' @rdname spregimes
#' @method print summary.spregimes
#' @import stats
#' @export
print.summary.spregimes <- function(x,
                                    digits = max(5, getOption("digits") - 3),
                                    ...)
{
  if("error_regimes" %in% class(x)){

    if(!is.null(unlist(x[[3]]))){
      if(is.null(x[[4]])){
        cat("\n ")
        cat("        --------------------------------------------\n")
        cat("                 Spatial Error Regimes Model      \n")
        cat("             and additional endogenous variables               \n")
        cat("        ---------------------------------------------\n")
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
        cat("\n ")
        cat("             --------------------------------------------\n")
        cat("                      Spatial Error Regimes Model \n")
        cat("                   with spatially lagged regressors \n")
        cat("                  and additional endogenous variables \n")
        cat("             --------------------------------------------\n")
        cat("\nCall:\n")
        cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

        cat("\nCoefficients:\n")
        printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)


        cat("\nEndogenous variables:\n")

        cat(paste(unlist(x[[3]]), sep=" "))

        cat("\nInstruments:\n")

        cat(paste(x[[4]], sep=" "))

      }
    }

    else{
      if(is.null(x[[4]])){
        cat("\n ")
        cat("               -----------------------------------\n")
        cat("                   Spatial Error Regimes Model       \n")
        cat("               -----------------------------------\n")
        cat("\nCall:\n")
        cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

        cat("\nCoefficients:\n")
        printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
      }
      else{
        cat("\n ")
        cat("                 ----------------------------------------\n")
        cat("                       Spatial Error Regimes Model        \n")
        cat("                     with spatially lagged regressors        \n")
        cat("                 ----------------------------------------\n")
        cat("\nCall:\n")
        cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

        cat("\nCoefficients:\n")
        printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

      }
    }

  }

  if("sarar_regimes" %in% class(x)){
    if(!is.null(x[[5]])){
      cat("\n ")
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
      cat("\n ")
      cat("        ------------------------------------------------------------\n")
      cat("                       Spatial SARAR Regimes Model       \n")
      cat("        ------------------------------------------------------------\n")
      cat("\nCall:\n")
      cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

      cat("\nCoefficients:\n")
      printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)


    }

  }

  if("lag_regimes" %in% class(x)){
    if(is.null(x[[5]])){

      if(is.null(x[[6]])){
        cat("\n ")
        cat("        ------------------------------------------------------------\n")
        cat("                          Spatial Lag Regimes Model \n")
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
        cat("\n ")
        cat("        ------------------------------------------------------------\n")
        cat("                          Spatial Durbin Regimes Model \n")
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

    }
    else{

      if(is.null(x[[6]])){
        cat("\n ")
        cat("        ------------------------------------------------------------\n")
        cat("                          Spatial Lag Regimes Model \n")
        cat("                     with additional endogenous variables \n")
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
        cat("\n ")
        cat("        ------------------------------------------------------------\n")
        cat("                          Spatial Durbin Regimes Model \n")
        cat("                     with additional endogenous variables \n")
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

    }

  }

  if("ols_regimes" %in% class(x)){


    if(!is.null(unlist(x[[3]]))){
      cat("\n")
      cat("        ------------------------------------------------------------\n")
      cat("                Regimes Model with spatially lagged regressors      \n")
      cat("                     and additional endogenous variables               \n")
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
      cat("\n")
      cat("        ------------------------------------------------------------\n")
      cat("                Regimes Model with spatially lagged regressors      \n")
      cat("        ------------------------------------------------------------\n")
      cat("\nCall:\n")
      cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

      cat("\nCoefficients:\n")
      printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)


    }

  }

  if("lm" %in% class(x)){
    cat("\n")
    cat("                 --------------------------------\n")
    cat("                           Regimes Model \n")
    cat("                 --------------------------------\n")

    cat("\nCall:\n")
    cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    cat("\nCoefficients:\n")
    printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

  }

  if("ivregimes" %in% class(x)){
    cat("\n")
    cat("                ------------------------------------\n")
    cat("                          IV Regimes Model \n")
    cat("                ------------------------------------\n")
    cat("\nCall:\n")
    cat(paste(deparse(x[[2]]), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    cat("\nCoefficients:\n")
    printCoefmat(x$CoefTable, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

    cat("\nEndogenous variables:\n")

    cat(paste(x[[3]], x[[4]], sep=" "))

    cat("\nInstruments:\n")

    cat(paste(x[[5]], sep=" "))
  }


   invisible(x)
}


#' @rdname spregimes
#' @method residuals spregimes
#' @export
residuals.spregimes <- function(object, ...){
  residuals <- object[[1]][[3]]
  if("lm" %in% class(object)) residuals <- residuals(object[[1]])
  return(residuals)
}

#' @rdname spregimes
#' @method fitted spregimes
#' @export
fitted.spregimes <- function(object, ...){
  fitted  <- object[[1]][[6]]
  if("lm" %in% class(object)) fitted <- fitted(object[[1]])
  return(fitted)
}




