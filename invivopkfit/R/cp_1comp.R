#'Analytical 1-compartment model
#'
#'Calculates plasma concentrations vs. time according to the analytical solution
#'for the 1-compartment model, for single bolus doses (IV and/or oral).
#'
#'
#'# Required parameters
#'
#'`params` must include the following named items:
#'   \describe{
#'   \item{kelim}{Elimination rate, 1/h.}
#'   \item{Vdist}{Apparent volume of central compartment, L/kg BW. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#'For oral administration (if any `iv.dose == FALSE`), `params` must also
#'include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_Vdist`}
#'   \item{kgutabs}{Rate of absorption from gut, 1/h.}
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
#'and `Fgutabs_V1`.)
#'
#'@param params A named list of numeric parameter values. See Details for
#'  requirements.
#'@param time A numeric vector of times in hours, reflecting the time point when
#'  concentration is measured after the corresponding single bolus dose. Must be
#'  same length as `dose` and `iv.dose`, or length 1.
#'@param dose A numeric vector of doses in mg/kg, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `iv.dose`, or
#'  length 1.
#'@param iv.dose A logical vector, reflecting the route of administration of
#'  each single bolus dose. TRUE for single IV bolus dose; FALSE for single oral
#'  bolus dose.  Must be same length as `time` and `dose`, or length 1.
#'
#'@return A vector of plasma concentration values (mg/L) corresponding to
#'  `time`.
#'
#'@author Caroline Ring, John Wambaugh
#'
#'@export cp_1comp
cp_1comp <- function(params, time, dose, iv.dose){

  #check whether lengths of time, dose, and iv.dose match
  time_len <- length(time)
  dose_len <- length(dose)
  ivdose_len <- length(iv.dose)
  len_all <- c(time_len, dose_len, ivdose_len)
  #Cases:
  # All three lengths are the same -- OK
  # Two lengths are the same and the third is 1 -- OK
  # Two lengths are 1 and the third is not 1 -- OK
  # Otherwise there is a problem
  good_len <- (length(unique(len_all)) == 1) |
    (length(unique(len_all)) == 2 &
       sum(len_all == 1) %in% c(1, 2))

  if(!good_len){
    stop(paste0("invivopkfit::cp_1comp(): ",
                "'time', 'dose', and 'iv.dose' ",
                "must either be the same length or length 1.\n",
                "'time' is length ", time_len, "\n",
                "'dose' is length ", dose_len, "\n",
                "'iv.dose' is length ", ivdose_len, "\n")
    )
  }

  #if any are length-1, repeat them to match the longest
  max_len <- max(len_all)
  time <- rep(time, length.out = max_len)
  dose <- rep(dose, length.out = max_len)
  iv.dose <- rep(iv.dose, length.out = max_len)

  #compute Fgutabs/Vdist if necessary
  if(all(c("Fgutabs", "Vdist") %in% names(params))){
    params$Fgutabs_Vdist <- params$Fgutabs/params$Vdist
  }

  #if Vdist and Fgutabs_Vdist are provided, but not Fgutabs, compute Fgutabs
  if(all(c("Vdist", "Fgutabs_Vdist") %in% names(params)) &
     !("Fgutabs" %in% names(params))
  ){
    params$Fgutabs <- params$Fgutabs_Vdist * params$Vdist
  }

  #drop any length-0 params
  param_length <- sapply(params, length)
  params <- params[param_length>0]

  if(any(iv.dose %in% FALSE)){

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
  }

  if(any(iv.dose %in% TRUE)){
    if(!all(c("kelim", "Vdist") %in% names(params))){
      stop(paste0("cp_1comp(): Error: For 1-compartment IV model, ",
                  "missing parameter(s): ",
                  paste(setdiff(c("kelim", "Vdist"), names(params)),
                        collapse = ", "),
      )
      )
    }
  }

  cp <- vector(mode = "numeric", length = length(time))

  #IV model\
  if(any(iv.dose %in% TRUE)){
  cp[iv.dose %in% TRUE] <- dose[iv.dose %in% TRUE]*
    exp(-params$kelim *
          time[iv.dose %in% TRUE])/
    params$Vdist
  }

  #Oral model
  if(any(iv.dose %in% FALSE)){

    if(params$kelim != params$kgutabs){
      #the usual case: kelim != kgutabs
      cp[iv.dose %in% FALSE] <- (params$Fgutabs_Vdist *
                                   dose[iv.dose %in% FALSE] *
                                   params$kgutabs)/
        (params$kgutabs - params$kelim) *
        (exp(-params$kelim * time[iv.dose %in% FALSE]) -
           exp(-params$kgutabs* time[iv.dose %in% FALSE]))

    }else{ #in case kelim = kgutabs, use the alternate model equation
      cp[iv.dose %in% FALSE] <- params$Fgutabs_Vdist *
        dose[iv.dose %in% FALSE] * params$kelim *
        time[iv.dose %in% FALSE] *
        exp(-params$kelim * time[iv.dose %in% FALSE])
    }
  }


  return(cp)
}
