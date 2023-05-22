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
#' - `param_name`: Names of the model parameters
#' - `param_units`: Units of the model parameters
#' - `optimize_param`: TRUE if each parameter is to be estimated from the data; FALSE otherrwise
#' - `use_param`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
#'@author Caroline Ring

get_params_flat <- function(data,
                            lower_bound = c(Vdist = 1e-4, #Vdist
                                            Fgutabs = 1e-4, #Fgutabs
                                            Fgutabs_Vdist = 1e-4, #Fgutabs_Vdist
                                            Rblood2plasma = 1e-4), #Rblood2plasma

                            upper_bound = c(Vdist = 1e6, #Vdist
                                            Fgutabs = 1, #Fgutabs
                                            Fgutabs_Vdist = 1e4, #Fgutabs_Vdist
                                            Rblood2plasma = 1e6),
                            param_units = ggplot2::aes(Vdist = paste0("(", #Vdist
                                                                      unique(Dose.Units),
                                                                      ")",
                                                                      "/",
                                                                      "(",
                                                                      unique(Conc.Units),
                                                                      ")"),
                                                       Fgutabs = "unitless fraction", #Fgutabs
                                                       Fgutabs_Vdist = paste0("(", #Fgutabs_Vdist
                                                                              unique(Conc.Units),
                                                                              ")",
                                                                              "/",
                                                                              "(",
                                                                              unique(Dose.Units),
                                                                              ")"),
                                                       Rblood2plasma = "unitless ratio")){
  #param names
  param_name <-c("Vdist",
                  "Fgutabs",
                  "Fgutabs_Vdist",
                  "Rblood2plasma")

  lower_bound_default = c(Vdist = 1e-4, #Vdist
                          Fgutabs = 1e-4, #Fgutabs
                          Fgutabs_Vdist = 1e-4, #Fgutabs_Vdist
                          Rblood2plasma = 1e-4) #Rblood2plasma

  lower_bound_missing <- setdiff(names(lower_bound_default),
                                 names(lower_bound))
  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]

  upper_bound_default = c(Vdist = 1e6, #Vdist
                          Fgutabs = 1, #Fgutabs
                          Fgutabs_Vdist = 1e4, #Fgutabs_Vdist
                          Rblood2plasma = 1e6)

  upper_bound_missing <- setdiff(names(upper_bound_default),
                                 names(upper_bound))
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]


  optimize_param <- rep(TRUE, length(param_name))

  use_param <- rep(TRUE, length(param_name))


if(!("oral" %in% data$Route)){
  #if no oral data, can't fit Fgutabs, or Fgutabs_Vdist,
  #and they won't be used.
  optimize_param[param_name %in% c("Fgutabs", "Fgutabs_Vdist")] <- FALSE
  use_param[param_name %in% c("Fgutabs", "Fgutabs_Vdist")] <- FALSE
}else{ #if yes oral data:
  if("iv" %in% data$Route){
    #if both oral and IV data, Fgutabs and Vdist can be fit separately, so turn off Fgutabs_Vdist
    optimize_param[param_name %in% c("Fgutabs_Vdist")] <- FALSE
    use_param[param_name %in% c("Fgutabs_Vdist")] <- FALSE
  }else{
    #if oral ONLY:
    #cannot fit Fgutabs and Vdist separately, so turn them off.
    optimize_param[param_name %in% c("Fgutabs", "Vdist")] <- FALSE
    use_param[param_name %in% c("Fgutabs", "Vdist")] <- FALSE
  }
}

  #if both "blood" and "plasma" are not in data,
  #then Rblood2plasma will not be optimized
  if(!(all(c("blood", "plasma") %in% data$Media))){
    optimize_param[param_name %in% "Rblood2plasma"] <- FALSE
  }

  # #if no medium is "blood" then Rblood2plasma will not be used at all
  # if(!("blood" %in% data$Media)){
  #   use_param[param_name %in% "Rblood2plasma"] <- FALSE
  #   #otherwise, if blood-only data, Rblood2plasma will be used, but held constant at 1
  # }

  #get param units based on data
  param_units_vect <- sapply(param_units,
                             function(x) rlang::eval_tidy(x,
                                                          data = data),
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  param_units_vect <- param_units_vect[param_name]


  par_DF <- data.frame("param_name" = param_name,
                       "param_units" = param_units_vect,
                       "optimize_param" = optimize_param,
                       "use_param" = use_param,
                       "lower_bound" = lower_bound,
                       "upper_bound" = upper_bound)

  rownames(par_DF) <- par_DF$param_name

  par_DF <- get_starts_flat(data = data,
                            par_DF = par_DF)

  #check to ensure starting values are within bounds
  #if not, then replace them by a value halfway between bounds
  start_low <- (par_DF$start < par_DF$lower_bound) %in% TRUE
  start_high <- (par_DF$start > par_DF$upper_bound) %in% TRUE
  start_nonfin <- !is.finite(par_DF$start)

  par_DF[start_low | start_high | start_nonfin,
         "start"] <- rowMeans(cbind(par_DF[start_low | start_high | start_nonfin,
                                           c("lower_bound",
                                             "upper_bound")]))


  return(par_DF)

}
