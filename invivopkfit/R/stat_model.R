#' Add PK model(s) to be fitted
#'
#' @param model The name(s) of models to be fitted. These should be the names of objects of class `pk_model`. Built-in options are `"flat"`,
#'   `"1comp"`, and `"2comp"`. You may add your own model by using [pk_model()].

stat_model.pk <- function(model = c("flat", "1comp", "2comp"),
                      ...){
  stat_model <- list(name = "model")
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  stat_model$value <- argg
  #set class
  class(stat_model) <- c(class(settings_data), "pkproto", "pk_stat")

  return(stat_model)
}

#'Create a new `pk_model` object
#'
#'# `conc_fun` requirements
#'
#'`conc_fun` should be a function that takes the following arguments, and
#'returns a numeric vector of predicted tissue concentrations:
#'
#' - `params`: A named list of parameter values
#' - `time`: A numeric vector of time values
#' - `dose`: A numeric vector of dose values. Currently, only a single bolus dose at time 0 is supported.
#' - `iv.dose`: A logical vector that is TRUE if the dose was administered intravenously; FALSE if the dose was administered orally
#' - `medium`: The tissue in which concentration is to be predicted. Currently, only `blood` and `plasma` are supported.
#'
#'See [cp_1comp()], [cp_2comp()], [cp_flat()] for examples.
#'
#'#`auc_fun` requirements
#'
#'`auc_fun` should be a function that takes the same arguments as `conc_fun`,
#'and returns a numeric vector of predicted tissue AUCs (area under the
#'concentration-time curve).
#'
#'See [auc_1comp()], [auc_2comp()], [auc_flat()] for examples.
#'
#'# `get_opts_fun` requirements
#'
#'`get_opts_fun` should be a function whose first argument is a `data.frame`,
#'which will be the pre-processed data using `invivopkfit` harmonized variable
#'names. It should then return a `data.frame` with the following variables:
#'
#' - `param_name`: Character vector, listing parameter names for the model
#' - `param_units`: Character vector, listing units of each model parameter
#' - `optimize_param`: Logical (TRUE/FALSE), whether each parameter is to be estimated given the available data
#' - `use_param`: Logical (TRUE/FALSE), whether each parameter is to be used in the model even if it is not estimated (i.e., if a parameter value is to be held constant while the others are estimated, then `optimize_param` should be FALSE but `use_param` should be TRUE)
#'
#'# `lower_fun` requirements
#'
#'`lower_fun` should be a function whose
#'
#'
#'@param name Character: The name of the model
#'@param params Character vector: Parameter names of the model
#'@param conc_fun Name of a function to predict tissue concentrations. See
#'  Details for requirements.
#'@param auc_fun Name of the function to predict AUC (area under the
#'  concentration-time curve). See Details for requirements.
#'@param params_fun Name of the function that produces a data frame listing all
#'  model parameters and whether each one is to be optimized and/or used, given
#'  the available data
#' @param lower_fun Name of the function that produces a named vector of lower bounds for each parameter
#' @param upper_fun Name of the function that produces a named vector of upper bounds for each parameter
#' @param start_fun Name of the function that produces a named vector of starting values for the optimizer for each parameter
#' @return An object of class `pk_model`. Effectively, a named list containing all of the arguments provided to this function.
#' @author Caroline Ring
#' @export
pk_model <- function(name,
                     params,
                     conc_fun,
                     auc_fun,
                     params_fun,
                     lower_fun,
                     upper_fun,
                     start_fun){
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_model$value <- argg
  #set class
  class(this_model) <- c(class(this_model), "pk_model")
  return(this_model)
}
