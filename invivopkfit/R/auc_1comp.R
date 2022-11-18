#' Analytic AUC for 1-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
#'
#' # Model equations
#'
#' ## IV dosing
#'
#' \deqn{\frac{D}{V_d k_e}  (1 - \exp(-t k_e))}
#'
#' ## Oral dosing
#'
#' \deqn{\frac{dose Fa ka}{V_d (k_a-k_e)} ((\frac{1}{k_e} - \frac{1}{k_a}) +
#' (\frac{-\exp(-t k_e)}{k_e} + \frac{\exp(-t k_a)}{k_a} )}
#'
#'@param time A vector of times in hours
#'@param params A named list of model parameter values. Must include:
#'\describe{
#'\item{kelim}{Elimination rate, 1/h}
#'\item{Vdist}{Volume of distribution, L/kg body weight}}
#'For oral administration (\code{iv.dose} FALSE), \code{params} must also include:
#'\describe{
#'\item{Fgutabs}{Oral bioavailability, unitless fraction}
#'\item{kgutabs}{Oral absorption rate, 1/h}}
#'@param dose Dose in mg/kg
#'@param iv.dose TRUE for single IV bolus dose; FALSE for single oral dose
#'
#'@return A vector of plasma AUC values, evaluated at each time point in `time`.
#'
#'@author Caroline Ring, John Wambaugh
#' @export auc_1comp

auc_1comp <- function(time,
                      params,
                      dose,
                      iv.dose){
  if (any(sapply(params,function(x) identical(x,numeric(0))))) return(0)

  if (is.null(params$Fgutabs)|is.na(params$Fgutabs))
  {
    params$Fgutabs<-1
  }
  if (is.null(params$kgutabs)|is.na(params$kgutabs))
  {
    params$kgutabs<-1
  }
  if (params$Fgutabs>1) params$Fgutabs<-1

  Vd <- params$Vdist
  ke <- params$kelim
  Fa <- params$Fgutabs
  ka <- params$kgutabs

  if (iv.dose){
    auc <- dose/(Vd*ke) * (1 - exp(-time*ke))
  }else{
    auc <- dose*Fa*ka/(Vd*(ka-ke)) *
      ((1/ke - 1/ka) +
      (-exp(-time*ke)/ke +
         exp(-time*ka)/ka
       ))
  }


  auc[auc<0] <- 0
  auc[!is.finite(auc)] <- 0

  return(auc)
}
