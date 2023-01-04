#'Time derivative of analytical 2-compartment model
#'
#'Calculates the time derivative (instantaneous rate of change) of plasma
#'concentration according to the analytical solution for the 2-compartment
#'model.
#'
#'This function is used by [postprocess_data()] to determine the time of peak
#'concentration for the 2-compartment model, by locating the point where the
#'time derivative of concentration crosses zero.
#'
#'@param params A named list of parameter values including the following:
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
#'  For oral administration (\code{iv.dose} FALSE), \code{params} must also
#'  include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   "Fgutabs_V1"}
#'   \item{kgutabs}{Rate of absorption from gut, 1/h.}
#'   }
#'
#'  For oral administration, in lieu of "V1" and "Fgutabs", you may instead
#'  provide "Fgutabs_V1", the ratio of Fgutabs to V1 (1/L). This is an alternate
#'  parameterization for situations where "Fgutabs" and "V1" are not
#'  identifiable separately (i.e. when oral data are available, but IV data are
#'  not). If "Fgutabs" and "V1" are provided, then "Fgutabs_V1" will not be
#'  used.
#'
#'@author Caroline Ring, John Wambaugh
#'@param time A numeric vector of times in hours, reflecting the time points
#'  when concentration is measured after the corresponding single bolus dose.
#'  Must be same length as `dose` and `iv.dose`, or length 1.
#'@param dose A numeric vector of doses in mg/kg, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `iv.dose`, or
#'  length 1.
#'@param iv.dose A logical vector, reflecting the route of administration of
#'  each single bolus dose. TRUE for single IV bolus dose; FALSE for single oral
#'  bolus dose. Must be same length as `time` and `dose`, or length 1.
#'@return A vector of instantaneous rates of change of plasma concentration
#'  values (mg/L/time) corresponding to each value in \code{time}
#'@export cp_2comp_dt
cp_2comp_dt <- function(params, time, dose, iv.dose)
{

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

  dcpdt <- vector(mode = "numeric", length = length(time))
  A <- vector(mode = "numeric", length = length(time))
  B <- vector(mode = "numeric", length = length(time))

  #see https://www.boomer.org/c/p4/c19/c1902.php
  #for these equations

  alpha_beta_sum <- params$kelim + params$k12 + params$k21
  alpha_beta_prod <- params$kelim * params$k21

  alpha <- (alpha_beta_sum + sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2
  beta <- (alpha_beta_sum - sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2


  A[iv.dose %in% TRUE] <- (dose[iv.dose %in% TRUE] *
                             (alpha - params$k21)) /
    (params$V1 * (alpha - beta))
  B[iv.dose %in% TRUE] <- (dose[iv.dose %in% TRUE] *
                             (params$k21 - beta)) /
    (params$V1 * (alpha - beta))

  dcpdt[iv.dose %in% TRUE] <-  A[iv.dose %in% TRUE] * - alpha *
    exp(-alpha * time[iv.dose %in% TRUE]) +
    B[iv.dose %in% TRUE] * - beta *
    exp(-beta * time[iv.dose %in% TRUE])


  #if any oral data, in which case params$kgutabs and params$Fgutabs_V1 exist:
  if(any(iv.dose %in% FALSE)){
  A[iv.dose %in% FALSE] <- (params$kgutabs * params$Fgutabs_V1 *
                              dose[iv.dose %in% FALSE] *
                              (alpha - params$k21)) /
    ( (params$kgutabs - alpha) * (alpha - beta))

  B[iv.dose %in% FALSE] <- (params$kgutabs * params$Fgutabs_V1 *
                              dose[iv.dose %in% FALSE] *
                              (params$k21 - beta)) /
    ( (params$kgutabs - beta) * (alpha - beta))

  dcpdt[iv.dose %in% FALSE] <-   A[iv.dose %in% FALSE] * - alpha *
    exp(-alpha * time[iv.dose %in% FALSE]) +
    B[iv.dose %in% FALSE] * - beta *
    exp(-beta * time[iv.dose %in% FALSE]) +
    -(A[iv.dose %in% FALSE] + B[iv.dose %in% FALSE]) *
    - params$kgutabs *
    exp(-params$kgutabs * time[iv.dose %in% FALSE])

  }

  return(dcpdt)
}
