#' Analytical AUC for the 2-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
#'
#' @param params A named list of parameter values including the following:
#'   \describe{
#'   \item{k12}{Rate at which compound moves from central to peripheral
#'   compartment, 1/h.}
#'   \item{k21}{Rate at which compound moves from peripheral to central
#'   compartment, 1/h.}
#'   \item{kelim}{Elimination rate, 1/h.}
#'   \item{V1}{Apparent volume of central compartment, L/kg BW. Or see below for
#'   "Fgutabs_V1"}
#'   }
#'
#'   For oral administration (\code{iv.dose} FALSE), \code{params} must also
#'   include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   "Fgutabs_V1"}
#'   \item{kgutabs}{Rate of absorption from gut, 1/h.}
#'   }
#'
#'   For oral administration, in lieu of "V1" and "Fgutabs", you may instead
#'   provide "Fgutabs_V1", the ratio of Fgutabs to V1 (1/L). This is an
#'   alternate parameterization for situations where "Fgutabs" and "V1" are not
#'   identifiable separately (i.e. when oral data are available, but IV data are
#'   not). If "Fgutabs" and "V1" are provided, then "Fgutabs_V1" will not be
#'   used.
#'
#' @param time A numeric vector of time values, in hours
#' @param dose A numeric vector of doses in mg/kg
#' @param iv.dose A logical vector: TRUE for single IV bolus dose, FALSE for single oral dose
#'@return A vector of plasma AUC values, evaluated at each time point in `time`.
#' @export auc_2comp
#' @author Caroline Ring, John Wambaugh
auc_2comp <- function(params, time, dose, iv.dose){

  if(all(c("Fgutabs", "V1") %in% names(params))){
    params$Fgutabs_V1 <- params$Fgutabs/params$V1
  }

  #drop any length-0 params
  param_length <- sapply(params, length)
  params <- params[param_length>0]

  if(any(iv.dose %in% FALSE)){
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

  if(any(iv.dose %in% TRUE)){
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
  }

  auc <- vector(mode = "numeric", length = length(time))
  A <- vector(mode = "numeric", length = length(time))
  B <- vector(mode = "numeric", length = length(time))

  #see https://www.boomer.org/c/p4/c19/c1902.php
  #for these equations

  alpha_beta_sum <- params$kelim + params$k12 + params$k21
  alpha_beta_prod <- params$kelim * params$k21

  alpha <- (alpha_beta_sum + sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2
  beta <- (alpha_beta_sum - sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2

  if(any(iv.dose %in% TRUE)){
  A[iv.dose %in% TRUE] <- (dose[iv.dose %in% TRUE] *
                             (alpha - params$k21)) /
    (params$V1 * (alpha - beta))
  B[iv.dose %in% TRUE] <- (dose[iv.dose %in% TRUE] *
                             (params$k21 - beta)) /
    (params$V1 * (alpha - beta))

auc[iv.dose %in% TRUE] <- A[iv.dose %in% TRUE]/alpha -
  A[iv.dose %in% TRUE]*exp(-time[iv.dose %in% TRUE]*alpha)/alpha +
  B[iv.dose %in% TRUE]/beta -
  B[iv.dose %in% TRUE]*exp(-time[iv.dose %in% TRUE]*beta)/beta
}

if(any(iv.dose %in% FALSE)){
  A[iv.dose %in% FALSE] <- (params$kgutabs * params$Fgutabs_V1 *
                              dose[iv.dose %in% FALSE] *
                              (alpha - params$k21)) /
    ( (params$kgutabs - alpha) * (alpha - beta))

  B[iv.dose %in% FALSE] <- (params$kgutabs * params$Fgutabs_V1 *
                              dose[iv.dose %in% FALSE] *
                              (params$k21 - beta)) /
    ( (params$kgutabs - beta) * (alpha - beta))

auc[iv.dose %in% FALSE] <- A[iv.dose %in% FALSE]/alpha -
  A[iv.dose %in% FALSE]*exp(-time[iv.dose %in% FALSE]*alpha)/alpha +
  B[iv.dose %in% FALSE]/beta -
  B[iv.dose %in% FALSE]*exp(-time[iv.dose %in% FALSE]*beta)/beta +
  (-A[iv.dose %in% FALSE] - B[iv.dose %in% FALSE])/params$kgutabs -
  (-A[iv.dose %in% FALSE] - B[iv.dose %in% FALSE])*
  exp(-time[iv.dose %in% FALSE]*params$kgutabs)/params$kgutabs
}

  return(auc)
}
