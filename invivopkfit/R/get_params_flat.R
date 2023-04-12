#' Get parameters to be optimized for flat model
#'
#' The full set of model parameters for the flat model includes
#' `Vdist`, and `Fgutabs`. Whether each one can be estimated
#' from the data depends on what routes of administration are included in the
#' data.
#'
#' ## IV data, no oral data
#'
#' If IV dosing data are available, but no oral dosing data are available, then
#' only the parameters `Vdist` and `kelim` will be estimated from the data. The
#' parameters `kgutabs` and `Fgutabs` cannot be estimated from IV data alone and will not be used to evaluate the model.
#'
#' ## Oral data, no IV data
#'
#' If oral dosing data are available, but no IV dosing data are available, then
#' the parameters `kelim` and `kgutabs` can be estimated from the data. However,
#' the parameters `Fgutabs` and `Vdist` cannot be identified separately. From
#' oral data alone, only the ratio `Fgutabs/Vdist` can be identified. This ratio
#' is represented by a single parameter named `Fgutabs_Vdist`. `Fgutabs` and
#' `Vdist` will not be estimated nor used in model evaluation, but `Fgutabs_Vdist` will be estimated, along
#' with `kelim` and `kgutabs`.
#'
#' ## Oral data and IV data
#'
#' If both oral and IV dosing data are available, then `Vdist`, `kelim`,
#' `kgutabs`, and `Fgutabs` will all be estimated from the data.
#'
#' # Blood and plasma data
#'
#' If both blood and plasma data are available, then `Rblood2plasma` will be estimated from the data.
#'
#' # Blood data, no plasma data
#'
#' If blood data but no plasma data are available, then `Rblood2plasma` will be used in model evaluation, but held constant at 1, not estimated from the data.
#'
#' # Plasma data, no blood data
#'
#' If plasma data but no blood data are available, then `Rblood2plasma` will neither be estimated nor be used in model evaluation at all.
#'
#'
#' @param data The data set to be fitted, with harmonized variable names (e.g. the result of [preprocess_data()]).
#'@return A `data.frame`with the following variables:
#' - `param_names`: Names of the model parameters
#' - `param_units`: Units of the model parameters
#' - `opt_params`: TRUE if each parameter is to be estimated from the data; FALSE otherrwise
#' - `use_params`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
#'@author Caroline Ring

get_params_flat <- function(data){
  #param names
  param_names <-c("Vdist",
                  "Fgutabs",
                  "Fgutabs_Vdist",
                  "Rblood2plasma")
 #param units
  param_units <- c(
                    paste0("(", #Vdist
                           unique(data$Dose.Units),
                           ")",
                           "/",
                           "(",
                           unique(data$Value.Units),
                           ")"),
                    "unitless fraction", #Fgutabs
                    paste0("(", #Fgutabs_Vdist
                           unique(data$Value.Units),
                           ")",
                           "/",
                           "(",
                           unique(data$Dose.Units),
                           ")"),
                    "unitless ratio" #Rblood2plasma
                   )

  opt_params <- rep(TRUE, length(param_names))

  use_params <- rep(TRUE, length(param_names))


if(!("po" %in% data$Route)){
  #if no oral data, can't fit Fgutabs, or Fgutabs_Vdist,
  #and they won't be used.
  opt_params[param_names %in% c("Fgutabs", "Fgutabs_Vdist")] <- FALSE
  use_params[param_names %in% c("Fgutabs", "Fgutabs_Vdist")] <- FALSE
}else{ #if yes oral data:
  if("iv" %in% data$Route){
    #if both oral and IV data, Fgutabs and Vdist can be fit separately, so turn off Fgutabs_Vdist
    opt_params[param_names %in% c("Fgutabs_Vdist")] <- FALSE
    use_params[param_names %in% c("Fgutabs_Vdist")] <- FALSE
  }else{
    #if oral ONLY:
    #cannot fit Fgutabs and Vdist separately, so turn them off.
    opt_params[param_names %in% c("Fgutabs", "Vdist")] <- FALSE
    use_params[param_names %in% c("Fgutabs", "Vdist")] <- FALSE
  }
}

  #if both "blood" and "plasma" are not in data,
  #then Rblood2plasma will not be optimized
  if(!(all(c("blood", "plasma") %in% data$Media))){
    opt_params[param_names %in% "Rblood2plasma"] <- FALSE
  }

  #if no medium is "blood" then Rblood2plasma will not be used at all
  if(!("blood" %in% data$Media)){
    use_params[param_names %in% "Rblood2plasma"] <- FALSE
    #otherwise, if blood-only data, Rblood2plasma will be used, but held constant at 1
  }

  par_DF <- data.frame("param_names" = param_names,
                       "param_units" = param_units,
                       "opt_params" = opt_params,
                       "use_params" = use_params)

  return(par_DF)

}
