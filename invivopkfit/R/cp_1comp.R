#'Analytical 1-compartment model
#'
#'Calculates plasma concentrations according to the analytical solution for the
#'1-compartment model.
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
#'@param time A vector of times in hours.
#'@param dose Dose in mg/kg
#'@param iv.dose Logical: TRUE for single IV bolus dose; FALSE for single oral
#'  dose
#'
#'@return A vector of plasma concentration values corresponding to `time`.
#'
#'@author Caroline Ring, John Wambaugh
#'
#'@export cp_1comp
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
