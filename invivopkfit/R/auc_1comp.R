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
                      iv.dose){

  #compute Fgutabs/Vdist if necessary
  if(all(c("Fgutabs", "Vdist") %in% names(params))){
    params$Fgutabs_Vdist <- params$Fgutabs/params$Vdist
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
        params$Fgutabs*params$kgutabs*
        (1/params$kgutabs - 1/params$kelim)/
        (params$Vdist*(-params$kelim + params$kgutabs)) +
        dose[iv.dose %in% FALSE]*
        params$Fgutabs*params$kgutabs*
        (exp(-time[iv.dose %in% FALSE]*params$kgutabs)/params$kgutabs -
           exp(-time[iv.dose %in% FALSE]*params$kelim)/params$kelim)/
        (params$Vdist*(-params$kelim + params$kgutabs))

    }else{ #in case kelim = kgutabs, use the alternate model equation
      auc[iv.dose %in% FALSE] <- dose[iv.dose %in% FALSE]*params$Fgutabs/
        (params$Vdist*params$kelim) +
        (-dose[iv.dose %in% FALSE]*params$Fgutabs*
           time[iv.dose %in% FALSE]*params$kelim -
           dose[iv.dose %in% FALSE]*params$Fgutabs)*
        exp(-time[iv.dose %in% FALSE]*params$kelim)/
        (params$Vdist*params$kelim)
    }
  }


  return(auc)
}
