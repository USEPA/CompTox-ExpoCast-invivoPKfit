#' Analytical 2-compartment model
#'
#' Calculates plasma concentration according to the analytical solution for the
#' 2-compartment model.
#'
#' @param params A named list of parameter values including the following:
#'   \describe{
#'   \item{k12}{Rate at which compound moves from central to peripheral
#'   compartment, 1/h.}
#'   \item{k21}{Rate at which compound moves from peripheral to central
#'   compartment, 1/h.}
#'   \item{kelim}{Elimination rate, 1/h.}
#'   \item{V1}{Apparent volume of central compartment, L. Or see below for
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
#' @author Caroline Ring, John Wambaugh
#' @param time A vector of time values, in hours
#' @param dose A dose in mg/kg
#' @param iv.dose TRUE for single IV bolus dose, FALSE for single oral dose
#' @return A vector of plasma concentration values corresponding to each value
#'   in \code{time}
#' @export cp_2comp
cp_2comp <- function(params, time, dose, iv.dose)
{

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

  cp <- A * exp(-alpha * time) + B * exp(-beta * time)
  }else{ #for oral dosing
  A <- (params$Fgutabs_V1 * dose * (alpha - params$k21)) / ( (alpha - beta))
  B <- (params$Fgutabs_V1 * dose * (params$k21 - beta)) / ((alpha - beta))
  C <- -(A + B)
  cp <-  A * exp(-alpha * time) + B * exp(-beta * time) + C * exp(-params$kgutabs * time)
  }

  return(cp)
}
