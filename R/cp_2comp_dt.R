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
#'  For oral administration (`route %in% "oral"`), \code{params} must also
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
#'  Must be same length as `dose` and `route`, or length 1.
#'@param dose A numeric vector of doses in mg/kg, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `route`, or length
#'  1.
#'@param route A character vector, reflecting the route of administration of
#'  each single bolus dose. Currently, only "iv" and "oral" are supported. Must
#'  be same length as `time` and `dose`, or length 1.
#'@param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as `time` and `dose`, or length 1.
#'@return A vector of instantaneous rates of change of plasma concentration
#'  values (mg/L/time) corresponding to each value in \code{time}
#'@export cp_2comp_dt
cp_2comp_dt <- function(params, time, dose, route, medium)
{

  #check whether lengths of time, dose, and route match
  time_len <- length(time)
  dose_len <- length(dose)
  route_len <- length(route)
  medium_len <- length(medium)

  len_all <- c(time_len, dose_len, route_len, medium_len)
  #Cases:
  # All three lengths are the same -- OK
  # Two lengths are the same and the third is 1 -- OK
  # Two lengths are 1 and the third is not 1 -- OK
  # Otherwise there is a problem
  good_len <- (length(unique(len_all)) == 1) |
    (length(unique(len_all)) == 2 &
       sum(len_all == 1) %in% c(1, 2))

  if(!good_len){
    stop(paste0("invivopkfit::cp_2comp_dt(): ",
                "'time', 'dose', 'route', and 'medium'",
                "must either be the same length or length 1.\n",
                "'time' is length ", time_len, "\n",
                "'dose' is length ", dose_len, "\n",
                "'route' is length ", route_len, "\n",
                "'medium' is length ", medium_len, "\n")
    )
  }

  #if any are length-1, repeat them to match the longest
  max_len <- max(len_all)
  time <- rep(time, length.out = max_len)
  dose <- rep(dose, length.out = max_len)
  route <- rep(route, length.out = max_len)
  medium <- rep(medium, length.out = max_len)

  if(all(c("Fgutabs", "V1") %in% names(params))){
    params$Fgutabs_V1 <- params$Fgutabs/params$V1
  }

  #if V1 and Fgutabs_V1 are provided, but not Fgutabs, compute Fgutabs
  if(all(c("V1", "Fgutabs_V1") %in% names(params)) &
     !("Fgutabs" %in% names(params))
  ){
    params$Fgutabs <- params$Fgutabs_V1 * params$V1
  }

  #if Fgutabs and Fgutabs_V1 provided, but not V1, compute V1
  if(all(c("Fgutabs", "Fgutabs_V1") %in% names(params)) &
     !("V1" %in% names(params))
  ){
    params$V1 <- params$Fgutabs / params$Fgutabs_V1
  }

  #drop any length-0 params
  param_length <- sapply(params, length)
  params <- params[param_length>0]

  if(any(route %in% "oral")){
    missing_params <- setdiff(c("kelim",
                                "k21",
                                "k12",
                                "Fgutabs_V1",
                                "kgutabs"),
                              names(params))
    if(length(missing_params)>0){
      stop(paste("cp_2comp_dt(): Error: For 2-compartment oral model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  if(any(route %in% "iv")){
    missing_params <- setdiff(c("kelim",
                                "V1",
                                "k12",
                                "k21"),
                              names(params))
    if(length(missing_params)>0){
      stop(paste("cp_2comp_dt(): Error: For 2-compartment IV model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  if(any(medium %in% "blood")){
    if(!("Rblood2plasma" %in% names(params))){
      stop(paste0("cp_2comp_dt(): Error: For 2-compartment model ",
                  "in blood: missing parameter Rblood2plasma"))
    }
  }

  dcpdt <- vector(mode = "numeric", length = max_len)
  A <- vector(mode = "numeric", length = max_len)
  B <- vector(mode = "numeric", length = max_len)

  #see https://www.boomer.org/c/p4/c19/c1902.php
  #for these equations

  alpha_beta_sum <- params$kelim + params$k12 + params$k21
  alpha_beta_prod <- params$kelim * params$k21

  alpha <- (alpha_beta_sum + sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2
  beta <- (alpha_beta_sum - sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2


  A[route %in% "iv"] <- (dose[route %in% "iv"] *
                             (alpha - params$k21)) /
    (params$V1 * (alpha - beta))
  B[route %in% "iv"] <- (dose[route %in% "iv"] *
                             (params$k21 - beta)) /
    (params$V1 * (alpha - beta))

  dcpdt[route %in% "iv"] <-  A[route %in% "iv"] * - alpha *
    exp(-alpha * time[route %in% "iv"]) +
    B[route %in% "iv"] * - beta *
    exp(-beta * time[route %in% "iv"])


  #if any oral data, in which case params$kgutabs and params$Fgutabs_V1 exist:
  if(any(route %in% "oral")){
  A[route %in% "oral"] <- (params$kgutabs * params$Fgutabs_V1 *
                              dose[route %in% "oral"] *
                              (alpha - params$k21)) /
    ( (params$kgutabs - alpha) * (alpha - beta))

  B[route %in% "oral"] <- (params$kgutabs * params$Fgutabs_V1 *
                              dose[route %in% "oral"] *
                              (params$k21 - beta)) /
    ( (params$kgutabs - beta) * (alpha - beta))

  dcpdt[route %in% "oral"] <-   A[route %in% "oral"] * - alpha *
    exp(-alpha * time[route %in% "oral"]) +
    B[route %in% "oral"] * - beta *
    exp(-beta * time[route %in% "oral"]) +
    -(A[route %in% "oral"] + B[route %in% "oral"]) *
    - params$kgutabs *
    exp(-params$kgutabs * time[route %in% "oral"])


  }

  if(any(medium %in% "blood")){
    dcpdt[medium %in% blood] <- params$Rblood2plasma * dcpdt[medium %in% "blood"]
  }

  return(dcpdt)
}
