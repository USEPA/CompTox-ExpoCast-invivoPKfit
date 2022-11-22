#'Analytical 1-compartment model
#'
#' Calculates plasma concentrations according to the analytical solution for the  1-compartment model.
#'
#'@param params A named list of model parameter values.
#'
#'Must include:
#'\describe{
#'\item{kelim}{Elimination rate, 1/h}
#'\item{Vdist}{Volume of distribution, L/kg body weight. Or see below for "Fgutabs_Vdist".}}
#'
#'For oral administration (\code{iv.dose} FALSE), \code{params} must also include:
#'\describe{
#'\item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for "Fgutabs_Vdist".}
#'\item{kgutabs}{Oral absorption rate, 1/h}
#'}
#'
#'For oral administration, in lieu of providing "Fgutabs" and "Vdist", you may
#'instead provide "Fgutabs_Vdist", the ratio of Fgutabs to Vdist (kg body
#'weight/L). This is an alternate parameterization for cases when Fgutabs and
#'Vdist are not separately identifiable, i.e. when only oral data are available
#'with no IV data. if "Fgutabs" and "Vdist" are provided, then "Fgutabs_Vdist"
#'will not be used.
#'}
#'@param time A vector of times in hours.
#'@param dose Dose in mg/kg
#'@param iv.dose Logical: TRUE for single IV bolus dose; FALSE for single oral dose
#'
#'@return A vector of plasma concentration values corresponding to `time`.
#'
#'@author Caroline Ring, John Wambaugh
#'
#' @export cp_1comp
cp_1comp <- function(params, time, dose, iv.dose){



  if (iv.dose){
    #for iv administration
    #check for needed params
    if(!all(c("kelim", "Vdist") %in% names(params))){
      stop(paste0("cp_1comp(): Error: For 1-compartment IV model, ",
      "missing parameter(s): ",
           paste(setdiff(c("kelim", "Vdist"), names(params)),
                 collapse = ", "),
      )
      )
    }

    cp <- dose*exp(-params$kelim * time)/params$Vdist

    }else{
      #for oral administration

      if(all(c("Fgutabs", "Vdist") %in% names(params))){
        params$Fgutabs_Vdist <- params$Fgutabs/params$Vdist
      }

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

      if(!(params$kelim == params$kgutabs)){
        #the usual case: kelim != kgutabs
      const <- (params$Fgutabs_Vdist * dose * params$kgutabs)/
        (params$kgutabs - params$kelim)
    cp <- const * (exp(-params$kelim * time) - exp(-params$kgutabs* time))

      }else{ #in case kelim = kgutabs, use the alternate model equation
        cp <- params$Fgutabs_Vdist * dose * params$kelim *
          time * exp(-params$kelim * time)
      }
    }

  return(cp)
}
