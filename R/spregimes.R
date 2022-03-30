##### Functions for spatial regimes ####
#' Estimation of spatial regime models
#' @name spregimes
#' @param formula a symbolic description of the model.
#' @param data the data of class \code{data.frame}.
#' @param listw a spatial weighting matrix
#' @param listw2 an (optional) spatial weighting matrix for sarar models
#' @param rgv variable to identify the regimes
#' @param initial.value initial value for the spatial error parameter
#' @param wy_rg default \code{wy_rg = TRUE}, the lagged dependent variable varies for regime
#' @param model   one of ("sarar", "lag", "error")
#' @param het heteroskedastic variance-covariance matrix
#' @param step1.c use with error and sarar
#' @param control list of controls for the minimization
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
#' @return An object of class ``\code{lag_regimes}'', or \code{sarar_regimes}, or \code{error_regimes}
#' @import Formula sphet stats Matrix
#' @export


spregimes <- function(formula, data = list(), listw, listw2 = NULL, rgv = NULL,
                      initial.value = 0.2, wy_rg = FALSE,
                      model = c("sarar", "lag", "error"), het = FALSE,
                      step1.c = FALSE, control = list()){




  cl = match.call()
  switch(match.arg(model),
         # sarar = sarargmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog,
         #                  instruments = instruments, lag.instr = lag.instr, initial.value = initial.value,
         #                  het = het, verbose = verbose, na.action = na.action,
         #                  step1.c = step1.c, control = control, HAC = HAC, cl = cl, Durbin = Durbin),
         lag = lag_regimes(formula = formula, data = data, listw = listw, rgv = rgv,
                      het = het, cl = cl, wy_rg = wy_rg),
         # #error = errorgmm(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog,
         #                  instruments = instruments, lag.instr = lag.instr, initial.value = initial.value,
         #                  het = het, verbose = verbose, na.action = na.action,
         #                  step1.c = step1.c, control = control, HAC = HAC, cl = cl, Durbin = Durbin),
         # #ivhac = laghac(formula = formula, data = data, listw = listw, listw2 = listw2, endog = endog,
         #                instruments = instruments, lag.instr = lag.instr,  verbose = verbose,
         #                na.action = na.action, het = het, HAC = HAC, distance = distance,
         #                type = type, bandwidth = bandwidth, cl = cl, Durbin = Durbin),
         stop("Argument model incorrectly specified")
  )

}







