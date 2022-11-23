#' Analytical AUC for the 2-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
#'
#' @param params A named list of parameter values including the following:
#'   \describe{ \item{k12}{Rate at which compound moves from central to
#'   peripheral compartment} \item{k21}{Rate at which compound moves from
#'   peripheral to central compartment} \item{kelim}{Elimination rate}
#'   \item{V1}{Apparent volume of central compartment} } For oral administration
#'   (\code{iv.dose} FALSE), \code{params} must also include: \describe{
#'   \item{Fgutabs}{Oral bioavailability} \item{kgutabs}{Rate of absorption from
#'   gut} }
#'
#' @author Caroline Ring, John Wambaugh
#' @param time A vector of time values, in hours
#' @param dose A dose in mg/kg
#' @param iv.dose TRUE for single IV bolus dose, FALSE for single oral dose
#'@return A vector of plasma AUC values, evaluated at each time point in `time`.
#' @export auc_2comp
auc_2comp <- function(params, time, dose, iv.dose){

  if(all(c("Fgutabs", "V1") %in% names(params))){
    params$Fgutabs_Vdist <- params$Fgutabs/params$Vdist
  }

  #Check for needed params
  if (iv.dose){
    missing_params <- setdiff(c("kelim",
                                "V1",
                                "k12",
                                "k21"),
                              names(params))
    if(length(missing_params)>0){
      stop(paste("cp_2comp(): Error: For 2-compartment IV model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }else{
    #check needed params for oral dose
    missing_params <- setdiff(c("kelim",
                                "k21",
                                "k12",
                                "Fgutabs_V1",
                                "kgutabs"),
                              names(params))
    if(length(missing_params)>0){
      stop(paste("cp_2comp(): Error: For 2-compartment oral model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  alpha_beta_sum <- params$kelim + params$k12 + params$k21
  alpha_beta_prod <- params$kelim * params$k21

  alpha <- (alpha_beta_sum + sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2
  beta <- (alpha_beta_sum - sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2


  if (iv.dose){ #for IV dosing
    A <- (dose * (alpha - params$k21)) / (params$V1 * (alpha - beta))
    B <- (dose * (params$k21 - beta))/(params$V1 * (alpha - beta))

    auc <- A/alpha * (1 - exp(-time*alpha)) +
      B/beta * (1 - exp(-time*beta))

  }else{ #for oral dosing
    A <- (params$Fgutabs_V1 * dose * (alpha - params$k21)) / ( (alpha - beta))
    B <- (params$Fgutabs_V1 * dose * (params$k21 - beta)) / ((alpha - beta))

    auc <- A/alpha * (1 - exp(-time*alpha)) +
      B/beta * (1 - exp(-time*beta)) +
      (-A - B)/params$kgutabs * (1 - exp(-time*params$kgutabs))
  }

  return(auc)
}
