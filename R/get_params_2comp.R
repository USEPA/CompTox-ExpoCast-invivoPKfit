#' Get parameters for 2-compartment model
#'
#' Get parameters for 2-compartment model and determine whether each is to be
#' estimated from the data
#'
#' The full set of model parameters for the 1-compartment model includes `V1`,
#' `kelim`, `k12`, `k21`, `kgutabs`, and `Fgutabs`.
#'
#' ## IV data, no oral data
#'
#' If IV dosing data are available, but no oral dosing data are available, then
#' only the parameters `V1`, `kelim`, `k12`, and `k21` will be estimated from
#' the data. The parameters `kgutabs` and `Fgutabs` cannot be estimated from IV
#' data alone.
#'
#' ## Oral data, no IV data
#'
#' If oral dosing data are available, but no IV dosing data are available, then
#' the parameters `kelim`, `k12`, `k21`, and `kgutabs` will be estimated from
#' the data. However, the parameters `Fgutabs` and `V1` cannot be identified
#' separately. From oral data alone, only the ratio `Fgutabs/V1` can be
#' identified. This ratio is represented by a single parameter named
#' `Fgutabs_V1`. `Fgutabs` and `V1` will not be optimized, but `Fgutabs_V1` will
#' be optimized, along with `kelim`, `k12`, `k21`, and `kgutabs`.
#'
#' ## Oral data and IV data
#'
#' If both oral and IV dosing data are available, then `V1`, `kelim`, `k12`,
#' `k21`, `kgutabs`, and `Fgutabs` will all be estimated from the data.
#'
#' # Blood and plasma data
#'
#' If both blood and plasma data are available, then `Rblood2plasma` will be estimated from the data.
#'
#'# Only one of blood or plasma data
#'
#'If only one of blood or plasma data are available, then `Rblood2plasma` will be
#'held constant at 1, not estimated from the data.
#'
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
#'@return A `data.frame`with the following variables:
#' - `param_name`: Character: Names of the model parameters
#' - `param_units`: Character: Units of the model parameters
#' - `optimize_param`: TRUE if each parameter is to be estimated from the data; FALSE otherwise
#' - `use_param`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
#' -`lower_bounds`: Numeric: The lower bounds for each parameter
#' - `upper_bounds`: Numeric: The upper bounds for each parameter
#' - `start`: Numeric: The starting guesses for each parameter
#'@author Caroline Ring
#' @family 2-compartment model functions
#' @family get_params functions
#' @family built-in model functions

get_params_2comp <- function(data,
                             lower_bound = ggplot2::aes(kelim = 0.5*log(2)/max(Time_trans),
                                                        k12 = 0.5*log(2)/max(Time_trans),
                                                        k21 = 0.5*log(2)/max(Time_trans),
                                                        V1 = 0.01,
                                                        Fgutabs = 0,
                                                        kgutabs = 0.5*log(2)/max(Time_trans),
                                                        Fgutabs_V1 = 0.01,
                                                        Rblood2plasma = 1e-2),
                             upper_bound = ggplot2::aes(kelim = 2*log(2)/min(Time_trans[Time_trans>0]),
                                                        k12 = 2*log(2)/min(Time_trans[Time_trans>0]),
                                                        k21 = 2*log(2)/min(Time_trans[Time_trans>0]),
                                                        V1 = 100,
                                                        Fgutabs = 1,
                                                        kgutabs = 2*log(2)/min(Time_trans[Time_trans>0]),
                                                        Fgutabs_V1 = 1e2,
                                                        Rblood2plasma = 100),
                             param_units = ggplot2::aes(kelim = paste0("1/", #kelim
                                                         unique(Time_trans.Units)),
                                          V1 = paste0("(", #V1
                                                         unique(Dose.Units),
                                                         ")",
                                                         "/",
                                                         "(",
                                                         unique(Conc.Units),
                                                         ")"),
                                          k21 = paste0("1/", #k21
                                                         unique(Time_trans.Units)),
                                          k12 = paste0("1/", #k12
                                                         unique(Time_trans.Units)),
                                          Fgutabs = "unitless fraction", #Fgutabs
                                          kgutabs = paste0("1/", #kgutabs
                                                           unique(Time_trans.Units)),
                                          Fgutabs_V1 = paste0("(", #Fgutabs_Vdist
                                                                 unique(Conc.Units),
                                                                 ")",
                                                                 "/",
                                                                 "(",
                                                                 unique(Dose.Units),
                                                                 ")"),
                                          Rblood2plasma = "unitless ratio")){
  #param names
  param_name <-c("kelim",
                  "V1",
                  "k21",
                  "k12",
                  "Fgutabs",
                  "kgutabs",
                  "Fgutabs_V1",
                  "Rblood2plasma")

  lower_bound_default = ggplot2::aes(kelim = 0.5*log(2)/max(Time_trans),
                                     k12 = 0.5*log(2)/max(Time_trans),
                                     k21 = 0.5*log(2)/max(Time_trans),
                                     V1 = 0.01,
                                     Fgutabs = 0,
                                     kgutabs = 0.5*log(2)/max(Time_trans),
                                     Fgutabs_V1 = 0.01,
                                     Rblood2plasma = 1e-2)

  lower_bound_missing <- setdiff(names(lower_bound_default),
                                 names(lower_bound))
  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]

  upper_bound_default =ggplot2::aes(kelim = 2*log(2)/min(Time_trans[Time_trans>0]),
                                    k12 = 2*log(2)/min(Time_trans[Time_trans>0]),
                                    k21 = 2*log(2)/min(Time_trans[Time_trans>0]),
                                    V1 = 100,
                                    Fgutabs = 1,
                                    kgutabs = 2*log(2)/min(Time_trans[Time_trans>0]),
                                    Fgutabs_V1 = 1e2,
                                    Rblood2plasma = 100)

  upper_bound_missing <- setdiff(names(upper_bound_default),
                                 names(upper_bound))
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]


  optimize_param <- rep(TRUE, length(param_name))

  use_param <- rep(TRUE, length(param_name))


if(!("oral" %in% data$Route)){
  #if no oral data, can't fit kgutabs, Fgutabs, or Fgutabs_V1,
  #and they won't be used.
  optimize_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_V1")] <- FALSE
  use_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_V1")] <- FALSE
}else{ #if yes oral data:
  if("iv" %in% data$Route){
    #if both oral and IV data, Fgutabs and V1 can be fit separately, so turn off Fgutabs_Vdist
    optimize_param[param_name %in% c("Fgutabs_V1")] <- FALSE
    use_param[param_name %in% c("Fgutabs_V1")] <- FALSE
  }else{
    #if oral ONLY:
    #cannot fit Fgutabs and V1 separately, so turn them off.
    optimize_param[param_name %in% c("Fgutabs", "V1")] <- FALSE
    use_param[param_name %in% c("Fgutabs", "V1")] <- FALSE
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
  par_DF <- get_starts_2comp(data = data,
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
