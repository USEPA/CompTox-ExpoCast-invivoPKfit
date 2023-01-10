#'Analytical 2-compartment model
#'
#'Calculates plasma concentration according to the analytical solution for the
#'2-compartment model.
#'
#'`params` must include the following named list elements:
#'   \describe{
#'   \item{k12}{Rate at which compound moves from central to peripheral
#'   compartment, 1/h.}
#'   \item{k21}{Rate at which compound moves from peripheral to central
#'   compartment, 1/h.}
#'   \item{kelim}{Elimination rate, 1/h.}
#'   \item{V1}{Apparent volume of central compartment, L/kg BW. Or see below for
#'   `Fgutabs_V1`}
#'   }
#'
#'For oral administration (any `iv.dose == FALSE`), `params` must also include
#'the following named list items:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_V1`}
#'   \item{kgutabs}{Rate of absorption from gut, 1/h.}
#'   }
#'
#'For oral administration, in lieu of `V1` and `Fgutabs`, you may instead
#'provide `Fgutabs_V1`, the ratio of Fgutabs to V1 (1/L). This is an alternate
#'parameterization for situations where `Fgutabs` and `V1` are not identifiable
#'separately (i.e. when oral TK data are available, but IV data are not). If
#'`Fgutabs` and `V1` are provided, they will override any value also provided
#'for `Fgutabs_V1`.
#'
#'If both oral and IV administration are specified (i.e., some `iv.dose == TRUE`
#'and some `iv.dose == FALSE`), then `V1` is required along with either
#'`Fgutabs` or `Fgutabs_V1`. (If `V1` and `Fgutabs_V1` are provided, but
#'`Fgutabs` is not provided, then `Fgutabs` will be calculated from `V1` and
#'`Fgutabs_V1`.)
#'
#'#'If `any(medium %in% 'blood')`, then `params` must also include
#'`Rblood2plasma`, the ratio of chemical concentration in whole blood to the
#'chemical concentration in blood plasma.
#'
#'@param params A named list of parameter values. See Details for requirements.
#'@param time A numeric vector of times in hours, reflecting the time points
#'  when concentration is measured after the corresponding single bolus dose.
#'  Must be same length as other arguments, or length 1.
#'@param dose A numeric vector of doses in mg/kg, reflecting single bolus doses
#'  administered at time 0. Must be same length as other arguments, or
#'  length 1.
#'@param iv.dose A logical vector, reflecting the route of administration of
#'  each single bolus dose. TRUE for single IV bolus dose; FALSE for single oral
#'  bolus dose. Must be same length as other arguments, or length 1.
#'@param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as other arguments, or length 1.
#'@return A vector of blood or plasma concentration values (mg/L) corresponding to each
#'  value in \code{time}
#'@export cp_2comp
#'
#'@author Caroline Ring, John Wambaugh
cp_2comp <- function(params, time, dose, iv.dose, medium = "plasma")
{

  #check whether lengths of time, dose, and iv.dose match
  time_len <- length(time)
  dose_len <- length(dose)
  ivdose_len <- length(iv.dose)
  medium_len <- length(medium)

  len_all <- c(time_len, dose_len, ivdose_len, medium_len)
  #Cases:
  # All three lengths are the same -- OK
  # Two lengths are the same and the third is 1 -- OK
  # Two lengths are 1 and the third is not 1 -- OK
  # Otherwise there is a problem
  good_len <- (length(unique(len_all)) == 1) |
    (length(unique(len_all)) == 2 &
       sum(len_all == 1) %in% c(1, 2))

  if(!good_len){
    stop(paste0("invivopkfit::cp_2comp(): ",
                "'time', 'dose', 'iv.dose', and 'medium'",
                "must either be the same length or length 1.\n",
                "'time' is length ", time_len, "\n",
                "'dose' is length ", dose_len, "\n",
                "'iv.dose' is length ", ivdose_len, "\n",
                "'medium' is length ", medium_len, "\n")
    )
  }

  #if any are length-1, repeat them to match the longest
  max_len <- max(len_all)
  time <- rep(time, length.out = max_len)
  dose <- rep(dose, length.out = max_len)
  iv.dose <- rep(iv.dose, length.out = max_len)
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

  #drop any length-0 params
  param_length <- sapply(params, length)
  params <- params[param_length>0]

  #check for any missing parameters
  #required params for oral dose
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

  #required params for IV dose
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

  if(any(medium %in% "blood")){
    if(!("Rblood2plasma" %in% names(params))){
      stop(paste0("cp_2comp(): Error: For 1-compartment model ",
                  "in blood: missing parameter Rblood2plasma"))
    }
  }

  cp <- vector(mode = "numeric", length = length(time))
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

  cp[iv.dose %in% TRUE] <-  A[iv.dose %in% TRUE] *
    exp(-alpha * time[iv.dose %in% TRUE]) +
    B[iv.dose %in% TRUE] *
    exp(-beta * time[iv.dose %in% TRUE])
  }


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

  cp[iv.dose %in% FALSE] <-   A[iv.dose %in% FALSE] *
    exp(-alpha * time[iv.dose %in% FALSE]) +
    B[iv.dose %in% FALSE] *
    exp(-beta * time[iv.dose %in% FALSE]) +
    -(A[iv.dose %in% FALSE] + B[iv.dose %in% FALSE]) *
    exp(-params$kgutabs * time[iv.dose %in% FALSE])

  }

  cp[medium %in% "blood"] <- params$Rblood2plasma * cp[medium %in% "blood"]

  return(cp)
}
