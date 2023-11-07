#' @title Estimation of spatial regimes models
#' @description The function \code{spregimes} deals
#' with the estimation of spatial regimes models.
#' This is a general function that allows the estimation
#' of various spatial specifications, including the spatial lag regimes model,
#' the spatial error regimes model, and the spatial SARAR regimes model.
#' Since the estimation is based on generalized method of moments (GMM),
#' endogenous variables can be included.
#' For further information on estimation, see details.
#'
#' @name spregimes
#' @param formula a symbolic description of the model of
#' the form \code{y ~ x_f | x_v | wx | h_f | h_v | wh}
#' where \code{y} is the dependent variable,
#' \code{x_f} are the regressors that do not vary by regimes,
#' \code{x_v} are the regressors that vary by regimes,
#' \code{wx} are the spatially lagged regressors,
#' \code{h_f} are the instruments that do not vary by regimes,
#' \code{h_v} are the instruments that vary by regimes,
#' \code{wh} are the spatially lagged instruments.
#' @param data the data of class \code{data.frame}.
#' @param model should be one of \code{c("sarar", "lag", "error", "ols")}
#' @param listw a spatial weighting matrix of class \code{listw}, \code{matrix} or \code{Matrix}
#' @param wy_rg default \code{wy_rg = FALSE}, the lagged dependent variable does not vary by regime (see details)
#' @param weps_rg default \code{weps_rg = FALSE}, if \code{TRUE} the spatial error term varies by regimes (see details)
#' @param initial.value initial value for the spatial error parameter
#' @param rgv an object of class \code{formula} to identify the regime variables
#' @param het heteroskedastic variance-covariance matrix
#' @param verbose print a trace of the optimization
#' @param control select arguments for the optimization
#' @param object an object of class spregimes
#' @param ... additional arguments
#' @param x an object of class spregimes
#' @param digits number of digits
#'
#' @examples
#' data("natreg")
#' data("ws_6")
#'
#' form <-  HR90  ~ 0 | MA90 + PS90 +
#' RD90 + UE90 | 0 | 0 | MA90 + PS90 +
#' RD90 + FH90 + FP89 + GI89 | 0
#'
#' form1 <-  HR90  ~ MA90 -1 |  PS90 +
#' RD90 + UE90 | 0 | MA90 -1 |  PS90 +
#' RD90 + FH90 + FP89 + GI89 | 0
#'
#' form2 <-  HR90  ~ MA90 -1 |  PS90 +
#' RD90 + UE90 | MA90 | MA90 -1 |  PS90 +
#' RD90 + FH90 + FP89 + GI89 | 0
#'
#' form3 <-  HR90  ~ MA90 -1 |  PS90 +
#' RD90 + UE90 | MA90 | MA90 -1 |  PS90 +
#' RD90 + FH90 + FP89 + GI89 | GI89
#'
#' form4 <-  HR90  ~ MA90 -1 |  PS90 +
#' RD90 + UE90 | MA90 + RD90 | MA90 -1 |  PS90 +
#' RD90 + FH90 + FP89 + GI89 | GI89
#'
#'
#' split  <- ~ REGIONS
#'
#' ###################################################
#' # Linear model with regimes and lagged regressors #
#' ###################################################

#' mod <- spregimes(formula = form2, data = natreg,
#' rgv = split, listw = ws_6, model = "ols")
#' summary(mod)
#'
#' mod1 <- spregimes(formula = form3, data = natreg,
#' rgv = split, listw = ws_6, model = "ols")
#' summary(mod1)
#'
#' mod2 <- spregimes(formula = form4, data = natreg,
#' rgv = split, listw = ws_6, model = "ols")
#' summary(mod2)
#'
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
#' Piras, G., Sarrias, M. (2023) Heterogeneous spatial models in R: spatial regimes models. J Spat Econometrics 4, 4. https://doi.org/10.1007/s43071-023-00034-1
#'
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
#'                                                  Gianfranco Piras, Paolo Postiglione (2022).  A deeper look at impacts in spatial Durbin model with sphet. \emph{Geographical Analysis}, 54(3), 664-684.
#'
#'                                                  Luc Anselin, Sergio J. Rey (2014). \emph{Modern Spatial Econometrics in Practice: A Guide to GeoDa, GeoDaSpace and PySal.} GeoDa Press LLC.
#'

#' @details
#' The function \code{spregimes} is a wrapper that allows the
#' estimation of a general
#' spatial regimes model.
#' For convenience and without loss of generality,
#' we assume the presence of only two regimes.
#' In this case the general model can be written as:
#' \deqn{
#' \begin{aligned}
#' y = &	W\begin{bmatrix}
#' y_1& 0 \\
#' 0 & y_2 \\
#' \end{bmatrix}
#' \begin{bmatrix}
#' \lambda_1 \\
#' \lambda_2 \\
#' \end{bmatrix}
#' +
#'  \begin{bmatrix}
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
#' + Y\pi +  \\
#' &
#' W\begin{bmatrix}
#' X_1& 0 \\
#' 0 & X_2 \\
#' \end{bmatrix}
#' \begin{bmatrix}
#' \delta_1 \\
#' \delta_2 \\
#' \end{bmatrix}+ WX\delta+
#'  W
#' \begin{bmatrix}
#' Y_1& 0 \\
#' 0 & Y_2 \\
#' \end{bmatrix}
#' \begin{bmatrix}
#' \theta_1 \\
#' \theta_2 \\
#' \end{bmatrix}
#' + WY\theta
#' +
#' \begin{bmatrix}
#' \varepsilon_1 \\
#' \varepsilon_2 \\
#' \end{bmatrix}
#' \end{aligned}
#' }
#' where
#' \deqn{
#' \begin{bmatrix}
#' \varepsilon_1 \\
#' \varepsilon_2 \\
#' \end{bmatrix}
#' =W	\begin{bmatrix}
#' \varepsilon_1&0 \\
#' 0&\varepsilon_2 \\
#' \end{bmatrix}
#' \begin{bmatrix}
#' \rho_1 \\
#' \rho_2 \\
#' \end{bmatrix}
#' +u \nonumber
#' }
#' The model includes the spatial lag of the dependent variable,
#' the spatial lag of the regressors,
#' the spatial lag of the errors
#' and, possibly, additional endogenous variables.
#' The function
#' \code{spregimes} estimates all of the nested
#' specifications deriving from this model.
#' There are, however, some restrictions.
#' For example, if \code{weps_rg} is set to TRUE,
#' all the regressors in the model should also vary by regimes.
#' The estimation of the different models relies heavily
#' on code available from the package \pkg{sphet}.
#'
#'\enumerate{
#' \item For the spatial lag (or Durbin) regimes model (i.e, when \eqn{\rho_1}
#' and \eqn{\rho_2} are zero), an instrumental variable
#' procedure is adopted, where the matrix of instruments
#' is formed by the spatial lags of the exogenous variables
#' and the additional instruments included in the \code{formula}.
#' A robust estimation
#' of the variance-covariance matrix can be obtained
#' by setting \code{het = TRUE}.
#'
#' \item For the spatial error regime models (i.e, when \eqn{\lambda_1}
#' and \eqn{\lambda_2} are zero), the spatial
#' coefficient(s)
#' are estimated with the GMM procedure described
#' in Kelejian and Prucha (2010) and Drukker et al., (2013).
#' The difference between Kelejian and Prucha (2010) and Drukker et al., (2013),
#' is that the former assume heteroskedastic innovations (\code{het = TRUE}),
#' while the latter does not (\code{het = FALSE}).
#'
#' \item For the SARAR regimes model, the estimation procedure
#' alternates a series of IV and GMM steps. The variance-covariance
#' can be estimated assuming that the innovations are homoskedastic (\code{het = FALSE})
#' as well as heteroskedastic (\code{het = TRUE}).
#' }
#'
#' @author Gianfranco Piras and Mauricio Sarrias
#' @return An object of class ``\code{spregimes}''
#' @import Formula sphet stats stringr
#' @rawNamespace import(Matrix,  except = c(cov2cor, toeplitz, update))
#' @export


spregimes <- function(formula, data = list(), model = c("sarar", "lag", "error", "ols"),
                      listw, wy_rg = FALSE, weps_rg = FALSE,
                       initial.value = NULL, rgv = NULL,
                      het = FALSE, verbose = FALSE, control = list()){




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


