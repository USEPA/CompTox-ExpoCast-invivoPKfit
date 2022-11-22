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

  Vd <- params$Vdist
  ke <- params$kelim
  Fa <- params$Fgutabs
  ka <- params$kgutabs

  if (iv.dose){
    #for iv administration
    #check for needed params
    if(!all(c("kelim", "Vdist") %in% names(params))){
      stop(paste0("cp_1comp(): Error: For 1-compartment IV model, ",
                  "missing parameter(s): ",
                  paste(setdiff(c("kelim", "Vdist"), names(params)),
                        collapse = ", "),
      )
      )
    }

    auc <- dose/(params$Vdist*params$kelim) * (1 - exp(-time*params$kelim))
  }else{
    #for oral administration

    if(all(c("Fgutabs", "Vdist") %in% names(params))){
      params$Fgutabs_Vdist <- params$Fgutabs/params$Vdist
    }

    #check for needed params
    if(!all(c("kelim", "kgutabs", "Fgutabs_Vdist") %in% names(params))){
      stop(paste0("cp_1comp(): Error: For 1-compartment oral model, ",
                  "missing parameter(s): ",
                  paste(setdiff(c("kelim", "kgutabs", "Fgutabs_Vdist"),
                                names(params)
                  ),
                  collapse = ", ")
      )
      )
    }

    if(!(params$kelim == params$kgutabs)){
      auc <- dose*params$Fgutabs_Vdist*params$kgutabs/
        (params$kgutabs-params$kelim) *
        ((1/params$kelim - 1/params$kgutabs) +
           (-exp(-time*params$kelim)/params$kelim +
              exp(-time*params$kgutabs)/params$kgutabs
           ))
    }else{
      auc <- dose*params$Fgutabs_Vdist/params$kelim +
        (-dose*params$Fgutabs_Vdist*time*params$kelim -
           dose*params$Fgutabs_Vdist)*
        exp(-time*params$kelim)/params$kelim
    }
  }

  return(auc)
}
