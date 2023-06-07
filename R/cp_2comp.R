#'Analytical 2-compartment model
#'
#'Calculates plasma concentration according to the analytical solution for the
#'2-compartment model.
#'
#'`params` must include the following named list elements:
#'   \describe{
#'   \item{k12}{Rate at which compound moves from central to peripheral
#'   compartment, 1/time.}
#'   \item{k21}{Rate at which compound moves from peripheral to central
#'   compartment, 1/time.}
#'   \item{kelim}{Elimination rate, 1/time.}
#'   \item{V1}{Apparent volume of central compartment, volume/unit BW. Or see below for
#'   `Fgutabs_V1`}
#'   }
#'
#'For oral administration (any `route %in% "oral"`), `params` must also include
#'the following named list items:
#'   \describe{
#'   \item{Fgutabs}{Oral bioavailability, unitless fraction. Or see below for
#'   `Fgutabs_V1`}
#'   \item{kgutabs}{Rate of absorption from gut, 1/time.}
#'   }
#'
#'For oral administration, in lieu of `V1` and `Fgutabs`, you may instead
#'provide `Fgutabs_V1`, the ratio of Fgutabs to V1 (1/volume). This is an alternate
#'parameterization for situations where `Fgutabs` and `V1` are not identifiable
#'separately (i.e. when oral TK data are available, but IV data are not). If
#'`Fgutabs` and `V1` are provided, they will override any value also provided
#'for `Fgutabs_V1`.
#'
#'If both oral and IV administration are specified (i.e., some `route %in% "iv"`
#'and some `route %in% "oral"`), then `V1` is required along with either
#'`Fgutabs` or `Fgutabs_V1`. (If `V1` and `Fgutabs_V1` are provided, but
#'`Fgutabs` is not provided, then `Fgutabs` will be calculated from `V1` and
#'`Fgutabs_V1`.)
#'
#'#'If `any(medium %in% 'blood')`, then `params` must also include
#'`Rblood2plasma`, the ratio of chemical concentration in whole blood to the
#'chemical concentration in blood plasma.
#'
#'@param params A named numeric vector of parameter values. See Details for requirements.
#'@param time A numeric vector of times, reflecting the time points
#'  when concentration is measured after the corresponding single bolus dose.
#'  Must be same length as other arguments, or length 1.
#'@param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as other arguments, or
#'  length 1.
#'@param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#'@param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as other arguments, or length 1.
#'@return A vector of blood or plasma concentration values (mass chemical/volume media) corresponding to each
#'  value in \code{time}
#'@export cp_2comp
#'
#'@author Caroline Ring, John Wambaugh
#' @family built-in model functions
#' @family 2-compartment model functions
#' @family model concentration functions
cp_2comp <- function(params, time, dose, route, medium = "plasma")
{
  #fill any missing parameters with NAs, and impute Fgutabs_v1 from Fgutabs and
  #V1 if necessary
  params <- fill_params_2comp(params)

#check that required params are present
  check_msg <- check_params_2comp(params,
                             route,
                             medium)
  if(!(check_msg %in% "Parameters OK")){
    stop(paste("cp_2comp():",
               check_msg))
  }

  #get transformed parameters for 2-comp model
  trans_params <- transformed_params_2comp(params)

  #for readability, assign params to variables inside this function
  for(x in names(params)){
    assign(x, unname(params[x]))
  }

  #for readability, assign transformed params to variables inside this function
  for(x in names(trans_params)){
    assign(x, unname(trans_params[x]))
  }

#get predicted concentration
  cp <- dose * ifelse(route %in% "iv",
                 A_iv_unit *
                   exp(-alpha * time) +
                   B_iv_unit *
                   exp(-beta * time),

                 A_oral_unit *
                   exp(-alpha * time) +
                   B_oral_unit *
                   exp(-beta * time) +
                   -(A_oral_unit + B_oral_unit) *
                   exp(-kgutabs * time)
  )

  cp <- ifelse(medium %in% "blood",
               Rblood2plasma * cp,
               cp)

  return(cp)
}
