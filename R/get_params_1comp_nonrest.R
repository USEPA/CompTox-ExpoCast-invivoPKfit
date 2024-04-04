#' Get parameters for 1-compartment model (non-restrictive)
#'
#' Get parameters for 1-compartment model that assumes restrictive clearance and
#' determine whether each is to be estimated from the data
#'
#'The full set of model parameters for the 1-compartment model includes `Vdist`,
#'`kelim`, `kgutabs`, `Fgutabs`, and `Rblood2plasma`. Whether each one can be
#'estimated from the data depends on what routes of administration are included
#'in the data. This model assumes restrictive clearance and will hold kelim constant.
#'
#'# IV data, no oral data
#'
#'If IV dosing data are available, but no oral dosing data are available, then
#'only the parameters `Vdist` will be estimated from the data. The
#'parameters `kgutabs` and `Fgutabs` cannot be estimated from IV data alone, and
#'will not be used in evaluating the model.
#'
#'# Oral data, no IV data
#'
#'If oral dosing data are available, but no IV dosing data are available, then
#'the parameter `kgutabs` can be estimated from the data. However,
#'the parameters `Fgutabs` and `Vdist` cannot be identified separately. From
#'oral data alone, only the ratio `Fgutabs/Vdist` can be identified. This ratio
#'is represented by a single parameter named `Fgutabs_Vdist`. `Fgutabs` and
#'`Vdist` will not be used to evaluate the model nor be estimated from data, but
#'`Fgutabs_Vdist` will be estimated from data, along with `kgutabs`.
#'
#'# Oral data and IV data
#'
#'If both oral and IV dosing data are available, then `Vdist`,
#'`kgutabs`, and `Fgutabs` will all be estimated from the data.
#'
#'# Blood and plasma data
#'
#'If both blood and plasma data are available, then `Rblood2plasma` will be
#'estimated from the data.
#'
#'# Only one of blood or plasma data
#'
#'If only one of blood or plasma data are available, then `Rblood2plasma` will be
#'held constant at 1, not estimated from the data.
#'
#'# Default lower and upper bounds for each parameter
#'
#'## Default lower and upper bounds for `kelim` and `kgutabs`
#'
#' Default bounds for time constant `kgutabs` are set based on
#' the time scale of the available data.
#'
#' Default bounds for `kelim` are set based on the calculation of constant `kelim`
#' from _httk_.
#'
#' The lower bounds are based on the assumption that absorption
#' are very slow compared to the time scale of the study. Specifically, the
#' lower bounds assume that absorption half-lives are twice as
#' long as the duration of the available study data, or `2*max(Time_trans)`.
#' Under this assumption, the corresponding  absorption time
#' constants would be `log(2)/(2*max(Time_trans))`. Therefore, the default lower
#' bounds for`kgutabs` are `log(2)/(2*max(Time_trans))`.
#'
#' Upper bounds are based on the opposite assumption: that
#' absorption are very fast compared to the time scale of the study.
#' Specifically, the upper bounds assume that the and absorption
#' half-lives are half as long as the time of the first observation after time
#' 0, or `0.5*min(Time_trans[Time_trans>0])`. Under this asumption, the
#' corresponding absorption time constants would be
#' `log(2)/(0.5*min(Time_trans[Time_trans>0]))`. Therefore, the default lower
#' bounds for `kgutabs` are
#' `log(2)/(0.5*min(Time_trans[Time_trans>0]))`.
#'
#' ## Default lower and upper bounds for `Vdist`
#'
#' By default, the lower bound for `Vdist` is 0.01, and the upper bound for
#' `Vdist` is 100. These values were chosen based on professional judgment.
#'
#' ## Default lower and upper bounds for `Fgutabs`
#'
#' By default, the lower bound for `Fgutabs` is 0.0, and the upper bound for
#' `Fgutabs` is 1. These are simply the bounds of the physically-meaningful
#' range for a fraction.
#'
#' ## Default lower and upper bounds for `Fgutabs_Vdist`
#'
#' By default, the lower bound for the ratio `Fgutabs_Vdist` is 0.01, and the
#' upper bound is 100. These values were chosen based on professional judgment.
#'
#' ## Default lower and upper bounds for `Rblood2plasma`
#'
#' By default, the lower bound for the blood:plasma partition coefficient
#' `Rblood2plasma` is 0.01, and the upper bound is 100. These values were chosen
#' based on professional judgment.
#'
#' # Starting values for each parameter
#'
#' Starting values for each parameter (starting guesses for the numerical
#' optimizer) are derived from the data using [get_starts_1comp()].
#'
#' If the starting values returned by [get_starts_1comp()] fall outside the
#' bounds for any parameter(s), then the starting value will be reset to a value
#' halfway between the lower and upper bounds for that parameter.
#'
#'@param data The data set to be fitted (e.g. the result of [preprocess_data()])
#'@param lower_bound A mapping specified using a call to [ggplot2::aes()],
#'  giving the lower bounds for each variable, as expressions which may include
#'  variables in `data`.
#' @param upper_bound A mapping specified using a call to [ggplot2::aes()],
#'  giving the upper bounds for each variable, as expressions which may include
#'  variables in `data`.
#'@param param_units A mapping specified using a call to [ggplot2::aes()],
#'  giving the units for each variable, as expressions which may include
#'  variables in `data`.
#' @param restrictive_clearance Logical. When NULL (default) it estimates `kelim`
#' but when set to TRUE or FALSE it uses a constant definition of `kelim` by
#' dividing total clearance by volume of distribution calculated by `httk`.
#'@return A `data.frame`with the following variables:
#' - `param_name`: Character: Names of the model parameters
#' - `param_units`: Character: Units of the model parameters
#' - `optimize_param`: TRUE if each parameter is to be estimated from the data; FALSE otherwise
#' - `use_param`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
#' -`lower_bounds`: Numeric: The lower bounds for each parameter
#' - `upper_bounds`: Numeric: The upper bounds for each parameter
#' - `start`: Numeric: The starting guesses for each parameter
#'@author Caroline Ring
#' @family 1-compartment model functions
#' @family get_params functions
#' @family built-in model functions

get_params_1comp_nonrest <- function(
    data,
    lower_bound = ggplot2::aes(kelim = 0,
                               Vdist = 0.01,
                               Fgutabs = 0.0,
                               kgutabs = log(2)/(2*max(Time_trans)),
                               Fgutabs_Vdist = 0.01,
                               Rblood2plasma = 1e-2),
    upper_bound = ggplot2::aes(kelim = 1E7,
                               Vdist = 100,
                               Fgutabs = 1,
                               kgutabs = log(2)/(0.5*min(Time_trans[Time_trans>0])),
                               Fgutabs_Vdist = 1e2,
                               Rblood2plasma = 100,),
    param_units = ggplot2::aes(kelim = paste0("1/", #kelim
                                              unique(Time_trans.Units)),
                               Vdist = paste0("(", #Vdist
                                              unique(Dose.Units),
                                              ")",
                                              "/",
                                              "(",
                                              unique(Conc.Units),
                                              ")"),
                               Fgutabs = "unitless fraction", #Fgutabs
                               kgutabs = paste0("1/", #kgutabs
                                                unique(Time_trans.Units)),
                               Fgutabs_Vdist = paste0("(", #Fgutabs_Vdist
                                                      unique(Conc.Units),
                                                      ")",
                                                      "/",
                                                      "(",
                                                      unique(Dose.Units),
                                                      ")"),
                               Rblood2plasma = "unitless ratio")){

  #param names
  param_name <- c("kelim",
                  "Vdist",
                  "Fgutabs",
                  "kgutabs",
                  "Fgutabs_Vdist",
                  "Rblood2plasma")

  #Default lower bounds, to be used in case the user specified a non-default
  #value for the `lower_bound` argument, but did not specify expressions for all
  #parameters. Any parameters not specified in the `lower_bound` argument will
  #take their default lower bounds defined here. This should be the same as the
  #default value for the `lower_bound` argument.
  lower_bound_default = ggplot2::aes(kelim = 0,
                                     Vdist = 0.01,
                                     Fgutabs = 0,
                                     kgutabs = log(2)/(2*max(Time_trans)),
                                     Fgutabs_Vdist = 0.01,
                                     Rblood2plasma = 1e-2)
  #which parameters did not have lower bounds specified in the `lower_bound`
  #argument?
  lower_bound_missing <- setdiff(names(lower_bound_default),
                                 names(lower_bound))
  #fill in the default lower bounds for any parameters that don't have them
  #defined in the `lower_bound` argument
  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]

  #Default upper bounds, to be used in case the user specified a non-default
  #value for the `upper_bound` argument, but did not specify expressions for all
  #parameters. Any parameters not specified in the `upper_bound` argument will
  #take their default upper bounds defined here. This should be the same as the
  #default value for the `upper_bound` argument.
  upper_bound_default = ggplot2::aes(kelim = 1E6,
                                     Vdist = 100,
                                     Fgutabs = 1,
                                     kgutabs = log(2)/(0.5*min(Time_trans[Time_trans>0])),
                                     Fgutabs_Vdist = 1e2,
                                     Rblood2plasma = 100)
  #which parameters did not have upper bounds specified in the `upper_bound`
  #argument?
  upper_bound_missing <- setdiff(names(upper_bound_default),
                                 names(upper_bound))
  #fill in the default upper bounds for any parameters that don't have them
  #defined in the `upper_bound` argument
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]

  #initialize optimization: start with optimize = TRUE for all params
  optimize_param <- rep(TRUE, length(param_name))

  #initialize whether each param is used: start with use = TRUE for all params
  use_param <- rep(TRUE, length(param_name))

  # Use constant kelim as calculated by httk
  optimize_param[param_name %in% "kelim"] <- FALSE

  #now follow the logic described in the documentation for this function:
  if(!("oral" %in% data$Route)){
    #if no oral data, can't fit kgutabs, Fgutabs, or Fgutabs_Vdist,
    #and they won't be used.
    optimize_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_Vdist")] <- FALSE
    use_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_Vdist")] <- FALSE
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

  param_units_vect <- sapply(param_units,
                             function(x) rlang::eval_tidy(x,
                                                          data = data),
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  param_units_vect <- param_units_vect[param_name]

  lower_bound_vect <- sapply(lower_bound,
                             function(x) rlang::eval_tidy(x,
                                                          data = data),
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  lower_bound_vect <- lower_bound_vect[param_name]

  upper_bound_vect <- sapply(upper_bound,
                             function(x) rlang::eval_tidy(x,
                                                          data = data),
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  upper_bound_vect <- upper_bound_vect[param_name]

  par_DF <- data.frame("param_name" = param_name,
                       "param_units" = param_units_vect,
                       "optimize_param" = optimize_param,
                       "use_param" = use_param,
                       "lower_bound" = lower_bound_vect,
                       "upper_bound" = upper_bound_vect)

  # now get starting values
  par_DF <-  get_starts_1comp_nonrest(data = data,
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
