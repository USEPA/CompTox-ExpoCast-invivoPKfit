#'Flat model
#'
#'Evaluates a "flat" model for concentration vs. time
#'
#'This function is used for model comparison: does a 1- or 2-compartment TK
#'model fit the data any better than this naive "flat" model?
#'
#'# Required parameters
#'
#'`params` must include the following named items:
#'   \describe{
#'   \item{Vdist}{Apparent volume of central compartment, L/kg BW. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#'For oral administration (if any `iv.dose == FALSE`), `params` must also
#'include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#'For oral administration, in lieu of `Vdist` and `Fgutabs`, you may instead
#'provide `Fgutabs_Vdist`, the ratio of Fgutabs to Vdist (1/L). This is an
#'alternate parameterization for situations where `Fgutabs` and `Vdist` are not
#'identifiable separately (i.e., when oral TK data are available, but IV data
#'are not). If `Fgutabs` and `Vdist` are provided, they will override any value
#'provided for `Fgutabs_Vdist`.
#'
#'If both oral and IV administration are specified (i.e., some `iv.dose == TRUE`
#'and some `iv.dose == FALSE`), then `Vdist` is required along with either
#'`Fgutabs` or `Fgutabs_Vdist`. (If `Vdist` and `Fgutabs_Vdist` are provided,
#'but `Fgutabs` is not provided, then `Fgutabs` will be calculated from `Vdist`
#'and `Fgutabs_Vdist`.)
#'
#'If `any(medium %in% 'blood')`, then `params` must also include
#'`Rblood2plasma`, the ratio of chemical concentration in whole blood to the
#'chemical concentration in blood plasma.
#'
#'@param params A named list of model parameter values. See Details for requirements.
#'@param time A numeric vector of times in hours.
#'@param dose A numeric vector of doses in mg/kg
#'@param iv.dose A logical vector: TRUE for single IV bolus dose; FALSE for
#'  single oral dose. Not used, but must be present for compatibility with other
#'  model functions.
#'@param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as other arguments, or length 1.
#'
#'@return A vector of plasma concentration values (mg/L) corresponding to
#'  \code{time}.
#'
#'@author Caroline Ring, John Wambaugh, Chris Cook
#'
#'@export cp_flat

cp_flat <- function(time, params, dose, iv.dose, medium) {

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
    stop(paste0("invivopkfit::cp_flat(): ",
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

  #check for required params

#if Fgutabs and Vdist provided, compute Fgutabs_Vdist
  if(all(c("Fgutabs", "Vdist")) %in% names(params)){
    params$Fgutabs_Vdist <- params$Fgutabs/params$Vdist
  }

  #if Vdist and Fgutabs_Vdist provided, but not Fgutabs, compute Fgutabs
  if(all(c("Vdist", "Fgutabs_Vdist") %in% names(params)) &
     !("Fgutabs" %in% names(params))
  ){
    params$Fgutabs <- params$Fgutabs_Vdist * params$Vdist
  }

  #if Fgutabs and Fgutabs_Vdist provided, but not Vdist, compute Vdist
  if(all(c("Fgutabs", "Fgutabs_Vdist") %in% names(params)) &
     !("Vdist" %in% names(params))
  ){
    params$Vdist <- params$Fgutabs / params$Fgutabs_Vdist
  }

  #drop any length-0 params
  param_length <- sapply(params, length)
  params <- params[param_length>0]

  if(any(iv.dose %in% TRUE)){
    missing_params <- setdiff(c("Vdist"),
                              names(params))
    if(length(missing_params)>0){
      stop(paste("cp_flat(): Error: For flat IV model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  if(any(iv.dose %in% FALSE)){
    missing_params <- setdiff(c("Fgutabs_Vdist"),
                              names(params))
    if(length(missing_params)>0){
      stop(paste("cp_flat(): Error: For flat PO model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  if(any(medium %in% "blood")){
    if(!("Rblood2plasma" %in% names(params))){
      stop(paste0("cp_flat(): Error: For flat model ",
                  "in blood: missing parameter Rblood2plasma"))
    }
  }

  cp <- vector(mode = "numeric", length = max_len)

  if(any(iv.dose %in% TRUE)){
  cp[iv.dose %in% TRUE] <- dose[iv.dose %in% TRUE]/params$Vdist
  }

  if(any(iv.dose %in% FALSE)){
  cp[iv.dose %in% FALSE] <- dose[iv.dose %in% FALSE]*params$Fgutabs_Vdist
  }

  if(any(medium %in% "blood")){
  cp[medium %in% "blood"] <- cp[medium %in% "blood"] * params$Rblood2plasma
  }

  return(cp)
}
