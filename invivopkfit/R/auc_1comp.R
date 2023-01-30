#' Analytic AUC for 1-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
#'
#'
#'@param params A named list of model parameter values. For either IV or oral
#'  dosing, this must include `kelim` (elimination rate, 1/h). For IV dosing, it
#'  must also include `Vdist`. For oral dosing, it must include `kgutabs` (oral
#'  absorption rate, 1/h), and either `Fgutabs` and `Vdist` (respectively,
#'  unitless fraction of oral dose that is absorbed, and volume of
#'  distribution), or `Fgutabs_Vdist` (the ratio of Fgutabs to Vdist).
#'  `Fgutabs_Vdist` is an alternate parameterization useful for parameter
#'  estimation when only oral data are available with no IV data, in which case
#'  only the ratio of `Fgutabs` to `Vdist` is identifiable. if `Fgutabs` and
#'  `Vdist` are provided along with `Fgutabs_Vdist`, then `Fgutabs_Vdist` will
#'  not be used.
#'@param time A numeric vector of times in hours.
#'@param dose A numeric vector of doses in mg/kg
#'@param iv.dose A logical vector: TRUE for single IV bolus dose; FALSE for single oral
#'  dose
#'
#'@return A vector of plasma AUC values (mg/L*time) corresponding to `time`.
#'
#'@author Caroline Ring, John Wambaugh
#' @export auc_1comp

auc_1comp <- function(params,
                      time,
                      dose,
                      iv.dose,
                      medium){

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
    stop(paste0("invivopkfit::auc_1comp(): ",
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

  #if Fgutabs and Fgutabs_Vdist provided, but not Vdist, compute Vdist
  if(all(c("Fgutabs", "Fgutabs_Vdist") %in% names(params)) &
     !("Vdist" %in% names(params))
  ){
    params$Vdist <- params$Fgutabs / params$Fgutabs_Vdist
  }

  #drop any length-0 params
  param_length <- sapply(params, length)
  params <- params[param_length>0]

  if(any(iv.dose %in% FALSE)){

    #check for needed params
    if(!all(c("kelim", "kgutabs", "Fgutabs_Vdist") %in% names(params))){
      stop(paste0("auc_1comp(): Error: For 1-compartment oral model, ",
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
      stop(paste0("auc_1comp(): Error: For 1-compartment IV model, ",
                  "missing parameter(s): ",
                  paste(setdiff(c("kelim", "Vdist"), names(params)),
                        collapse = ", ")
      )
      )
    }
  }

  if(any(medium %in% "blood")){
    if(!("Rblood2plasma" %in% names(params))){
      stop(paste0("auc_1comp(): Error: For 1-compartment model ",
                  "in blood: missing parameter Rblood2plasma"))
    }
  }

  auc <- vector(mode = "numeric", length = length(time))

  #IV model\
  if(any(iv.dose %in% TRUE)){
    auc[iv.dose %in% TRUE] <- dose[iv.dose %in% TRUE]/
      (params$Vdist*params$kelim) -
      dose[iv.dose %in% TRUE]*
      exp(-time[iv.dose %in% TRUE]*params$kelim)/
      (params$Vdist*params$kelim)
  }

  #Oral model
  if(any(iv.dose %in% FALSE)){

    if(params$kelim != params$kgutabs){
      #the usual case: kelim != kgutabs
      auc[iv.dose %in% FALSE] <- -dose[iv.dose %in% FALSE]*
        params$Fgutabs_Vdist*params$kgutabs*
        (1/params$kgutabs - 1/params$kelim)/
        ((-params$kelim + params$kgutabs)) +
        dose[iv.dose %in% FALSE]*
        params$Fgutabs_Vdist*params$kgutabs*
        (exp(-time[iv.dose %in% FALSE]*params$kgutabs)/params$kgutabs -
           exp(-time[iv.dose %in% FALSE]*params$kelim)/params$kelim)/
        ((-params$kelim + params$kgutabs))

    }else{ #in case kelim = kgutabs, use the alternate model equation
      auc[iv.dose %in% FALSE] <- dose[iv.dose %in% FALSE]*params$Fgutabs_Vdist/
        (params$kelim) +
        (-dose[iv.dose %in% FALSE]*params$Fgutabs_Vdist*
           time[iv.dose %in% FALSE]*params$kelim -
           dose[iv.dose %in% FALSE]*params$Fgutabs)*
        exp(-time[iv.dose %in% FALSE]*params$kelim)/
        (params$kelim)
    }
  }

  #convert blood concentrations using Rblood2plasma
  if(any(medium %in% "blood")){
    auc[medium %in% "blood"] <- params$Rblood2plasma * auc[medium %in% "blood"]
  }

  return(auc)
}
