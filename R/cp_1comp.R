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
#'   \item{kelim}{Elimination rate, 1/time.}
#'   \item{Vdist}{Apparent volume of central compartment, volume/unit BW. Or see below for
#'   `Fgutabs_Vdist`}
#'   }
#'
#'For oral administration (if any `route %in% "oral"`), `params` must also
#'include:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_Vdist`}
#'   \item{kgutabs}{Rate of absorption from gut, 1/time.}
#'   }
#'
#'For oral administration, in lieu of `Vdist` and `Fgutabs`, you may instead
#'provide `Fgutabs_Vdist`, the ratio of Fgutabs to Vdist (1/volume). This is an
#'alternate parameterization for situations where `Fgutabs` and `Vdist` are not
#'identifiable separately (i.e., when oral TK data are available, but IV data
#'are not). If `Fgutabs` and `Vdist` are provided, they will override any value
#'provided for `Fgutabs_Vdist`.
#'
#'If both oral and IV administration are specified (i.e., some `route %in% "iv"`
#'and some `route %in% "oral"`), then `Vdist` is required along with either
#'`Fgutabs` or `Fgutabs_Vdist`. (If `Vdist` and `Fgutabs_Vdist` are provided,
#'but `Fgutabs` is not provided, then `Fgutabs` will be calculated from `Vdist`
#'and `Fgutabs_Vdist`.)
#'
#'If `any(medium %in% 'blood')`, then `params` must also include
#'`Rblood2plasma`, the ratio of chemical concentration in whole blood to the
#'chemical concentration in blood plasma.
#'
#'@param params A named list of numeric parameter values. See Details for
#'  requirements.
#'@param time A numeric vector of times, reflecting the time point when
#'  concentration is measured after the corresponding single bolus dose. Must be
#'  same length as `dose` and `iv.dose`, or length 1.
#'@param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `iv.dose`, or
#'  length 1.
#'@param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#'@param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as `time` and `dose`, or length 1.
#'
#'@return A vector of blood or plasma concentration values  corresponding
#'  to `time`.
#'
#'@author Caroline Ring, John Wambaugh
#'
#'@export cp_1comp
#' @family built-in model functions
#' @family 1-compartment model functions
#' @family model concentration functions
cp_1comp <- function(params, time, dose, route, medium = 'plasma'){

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
    stop(paste0("invivopkfit::cp_1comp(): ",
                "'time', 'dose', 'iv.dose', and 'medium'",
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

  if(any(route %in% "oral")){

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

  if(any(route %in% "iv")){
    if(!all(c("kelim", "Vdist") %in% names(params))){
      stop(paste0("cp_1comp(): Error: For 1-compartment IV model, ",
                  "missing parameter(s): ",
                  paste(setdiff(c("kelim", "Vdist"), names(params)),
                        collapse = ", ")
      )
      )
    }
  }

  if(any(medium %in% "blood")){
    if(!("Rblood2plasma" %in% names(params))){
      stop(paste0("cp_1comp(): Error: For 1-compartment model ",
                  "in blood: missing parameter Rblood2plasma"))
    }
  }

  cp <- vector(mode = "numeric", length = max_len)

  #IV model\
  if(any(route %in% "iv")){
  cp[route %in% "iv"] <- dose[route %in% "iv"]*
    exp(-params$kelim *
          time[route %in% "iv"])/
    params$Vdist
  }

  #Oral model
  if(any(route %in% "oral")){

    if(params$kelim != params$kgutabs){
      #the usual case: kelim != kgutabs
      cp[route %in% "oral"] <- (params$Fgutabs_Vdist *
                                   dose[route %in% "oral"] *
                                   params$kgutabs)/
        (params$kgutabs - params$kelim) *
        (exp(-params$kelim * time[route %in% "oral"]) -
           exp(-params$kgutabs* time[route %in% "oral"]))

    }else{ #in case kelim = kgutabs, use the alternate model equation
      cp[route %in% "oral"] <- params$Fgutabs_Vdist *
        dose[route %in% "oral"] * params$kelim *
        time[route %in% "oral"] *
        exp(-params$kelim * time[route %in% "oral"])
    }
  }

#convert blood concentrations using Rblood2plasma
  if(any(medium %in% "blood")){
  cp[medium %in% "blood"] <- params$Rblood2plasma * cp[medium %in% "blood"]
  }

  return(cp)
}
