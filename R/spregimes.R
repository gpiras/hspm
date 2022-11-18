##### Functions for spatial regimes ####
#' @title Estimation of spatial regimes models
#' @name spregimes
#' @param formula a symbolic description of the model of the form \code{y ~ x_f | x_v | wx | h_f | h_v | wh}; where \code{y} is the dependent variable, \code{x_f} are the regressors that do not vary by regimes,  \code{x_v} are the regressors that vary by regimes, \code{wx} are the spatially lagged regressors, \code{h_f} are the instruments that do not vary by regimes,  \code{h_v} are the instruments that vary by regimes, \code{wh} are the spatially lagged instruments.
#' @param data the data of class \code{data.frame}.
#' @param listw a spatial weighting matrix of class \code{listw}, \code{matrix} or \code{Matrix}
#' @param rgv an object of class \code{formula} to identify the regime variables
#' @param initial.value initial value for the spatial error parameter
#' @param wy_rg default \code{wy_rg = FALSE}, the lagged dependent variable does not vary by regime (see details)
#' @param weps_rg default \code{weps_rg = FALSE}, the errors do not vary by regime (see details)
#' @param model one of \code{model = c("sarar", "lag", "error", "ols")}
#' @param het heteroskedastic variance-covariance matrix
#' @param control list of controls for the minimization
#' @param verbose print a trace of the optimization
#'
#' @examples
#' data("natreg")
#' data("ws_6")
#' form <-  HR90  ~ 0 | MA90 + PS90 +
#' RD90 + UE90 | 0 | 0 | MA90 + PS90 +
#' RD90 + FH90 + FP89 + GI89 | 0
#'
#' form1 <-  HR90  ~ MA90 -1 |  PS90 +
#' RD90 + UE90 | 0 | MA90 -1 |  PS90 +
#' RD90 + FH90 + FP89 + GI89 | 0
#'
#' split  <- ~ REGIONS
#'
#' ###############################
#' # Spatial Error regimes model #
#' ###############################
#' mod <- spregimes(formula = form, data = natreg,
#' rgv = split, listw = ws_6, model = "error", het = TRUE)
#' summary(mod)
#' mod1 <- spregimes(formula = form, data = natreg,
#' rgv = split, listw = ws_6, model = "error",
#' weps_rg = TRUE, het = TRUE)
#' summary(mod1)
#' mod2 <- spregimes(formula = form1, data = natreg,
#' rgv = split, listw = ws_6, model = "error", het = TRUE)
#' summary(mod2)
#'
#' ###############################
#' #  Spatial Lag regimes model  #
#' ###############################
#' mod4 <- spregimes(formula = form, data = natreg,
#' rgv = split, listw = ws_6, model = "lag",
#' het = TRUE, wy_rg = TRUE)
#' summary(mod4)
#' mod5 <- spregimes(formula = form1, data = natreg,
#' rgv = split, listw = ws_6, model = "lag",
#' het = TRUE, wy_rg = TRUE)
#' summary(mod5)
#'
#' ###############################
#' # Spatial SARAR regimes model #
#' ###############################
#' mod6 <- spregimes(formula = form, data = natreg,
#' rgv = split, listw = ws_6, model = "sarar",
#' het = TRUE, wy_rg = TRUE, weps_rg = TRUE)
#' summary(mod6)
#' mod7 <- spregimes(formula = form, data = natreg,
#' rgv = split, listw = ws_6, model = "sarar",
#' het = TRUE, wy_rg = FALSE, weps_rg = FALSE)
#' summary(mod7)
#' mod8 <- spregimes(formula = form1, data = natreg,
#' rgv = split, listw = ws_6, model = "sarar",
#' het = TRUE, wy_rg = TRUE, weps_rg = FALSE)
#' summary(mod8)
#' @references
#' Arraiz, I. and Drukker, M.D. and Kelejian, H.H. and Prucha, I.R. (2010)
#'   A spatial Cliff-Ord-type Model with Heteroskedastic Innovations: Small and Large Sample Results,
#'     \emph{Journal of Regional Sciences}, \bold{50}, pages 592--614.
#'
#'         Drukker, D.M. and Egger, P. and Prucha, I.R. (2013)
#'           On Two-step Estimation of a Spatial Auto regressive Model with Autoregressive
#'             Disturbances and Endogenous Regressors,
#'               \emph{Econometric Review}, \bold{32}, pages 686--733.
#'
#'
#'                 Kelejian, H.H. and Prucha, I.R. (2010)
#'                   Specification and Estimation of Spatial Autoregressive Models with Autoregressive and Heteroskedastic Disturbances,
#'                     \emph{Journal of Econometrics}, \bold{157}, pages 53--67.
#'
#'                                            Gianfranco Piras (2010). sphet: Spatial Models with Heteroskedastic Innovations in R. \emph{Journal of Statistical Software}, 35(1), 1-21. \doi{10.18637/jss.v035.i01}.
#'
#'                                                Roger Bivand, Gianfranco Piras (2015). Comparing Implementations of Estimation Methods for Spatial Econometrics. \emph{Journal of Statistical Software}, 63(18), 1-36. \doi{10.18637/jss.v063.i18}.
#'
#'                                                  Gianfranco Piras, Paolo Postiglione (2022).  A deeper look at impacts in spatial Durbin model with sphet. \emph{Geographical Analysis}, 54(3), 664-684. \url{https://onlinelibrary.wiley.com/doi/10.1111/gean.12318}
#'
#'                                                  Luc Anselin, Sergio J. Rey (2014). \emph{Modern Spatial Econometrics in Practice: A Guide to GeoDa, GeoDaSpace and PySal.} GeoDa Press LLC.
#'

#' @details
#' The general model contains the spatial lag of the dependent variable, the spatial lag of the regressors, the spatial lag of the errors and possibly additional endogenous variables.
#' \code{spregimes} estimate all of the nested specifications of this general model.
#' The regressors can be either "fixed" or varying by regime. However, if \code{weps_rg} is set to TRUE, all the regressors should vary by regime.
#'
#' @author Gianfranco Piras and Mauricio Sarrias
#' @return An object of class ``\code{lag_regimes}'', or \code{sarar_regimes}, or \code{error_regimes}
#' @import Formula sphet stats
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update))
#' @export


spregimes <- function(formula, data = list(), listw, rgv = NULL,
                      initial.value = NULL, verbose = FALSE, wy_rg = FALSE, weps_rg = FALSE,
                      model = c("sarar", "lag", "error", "ols"), het = FALSE,
                      control = list()){




  cl = match.call()
  switch(match.arg(model),
          sarar = sarar_regimes(formula = formula, data = data, listw = listw,
                           rgv = rgv, het = het, initial.value = initial.value,
                           verbose = verbose, control = control, cl = cl,
                           wy_rg = wy_rg, weps_rg = weps_rg),
         lag = lag_regimes(formula = formula, data = data, listw = listw, rgv = rgv,
                      het = het, cl = cl, wy_rg = wy_rg),
          error = error_regimes(formula = formula, data = data, listw = listw,  rgv = rgv,
                           initial.value = initial.value, het = het,
                           control = control, cl = cl, weps_rg = weps_rg, verbose = verbose),
         ols = ols_regimes(formula = formula, data = data, listw = listw, rgv = rgv,
                           het = het, cl = cl),
         stop("Argument model incorrectly specified")
  )

}








