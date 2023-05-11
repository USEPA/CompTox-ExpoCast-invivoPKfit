#' Add PK model(s) to be fitted
#'
#' @param model The name(s) of models to be fitted. These should be the names of objects of class `pk_model`. Built-in options are `"flat"`,
#'   `"1comp"`, and `"2comp"`. You may add your own model by using [pk_model()].

stat_model <- function(model = c("flat", "1comp", "2comp"),
                      ...){
  this_stat_model <- list()
  #pull pk_model objects by each name
for(this_model in model){
  this_stat_model[[this_model]] <- list()
  #check whether an object exists by this name
  if(exists(this_model)){
    this_model_obj <- get(this_model)
    #Check whether this is an object of class `pk_model`
    if(inherits(this_model_obj, "pk_model")){
      #if so, add it to the list
      this_stat_model[[this_model]] <- this_model_obj
    }
  }
}

  #set class
  class(this_stat_model) <- c(class(this_stat_model), "pkproto", "pk_stat_model")

  return(this_stat_model)
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
#' - `route`: A character vector of the route of administration. Currently, only `'oral'` and `'iv'` are supported.
#' - `medium`: The tissue in which concentration is to be predicted. Currently, only `'blood'` and `'plasma'` are supported.
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
#'# `params_fun` requirements
#'
#'`params_fun` should be a function whose first argument is a `data.frame`,
#'which will be the pre-processed data using `invivopkfit` harmonized variable
#'names. It may take additional arguments, which can be provided in
#'`params_fun_args`. The function should return a `data.frame` with the
#'following variables:
#'
#' - `param_name`: Character vector, listing parameter names for the model
#' - `param_units`: Character vector, listing units of each model parameter
#' - `optimize_param`: Logical (TRUE/FALSE), whether each parameter is to be estimated given the available data
#' - `use_param`: Logical (TRUE/FALSE), whether each parameter is to be used in the model even if it is not estimated (i.e., if a parameter value is to be held constant while the others are estimated, then `optimize_param` should be FALSE but `use_param` should be TRUE)
#' -`lower_bound`: Numerical. Lower bounds for each parameter. May be `-Inf` if no lower bound.  If `optimize_param` or `use_param` is FALSE, thenthe corresponding `lower_bound` will be ignored (because the parameter is not being estimated from the data).
#' - `upper_bound`: Numerical. Upper bounds for each parameter. May be `Inf` if no upper bound. If `optimize_param` or `use_param` is FALSE, then the corresponding `upper_bound` will be ignored (because the parameter is not being estimated from the data).
#' - `start`: Numerical. Starting values for estimating each parameter. If `optimize_param` is FALSE and `use_param` is TRUE, then the parameter will be held constant at the corresponding value in `start`. If `use_param` is FALSE, then the corresponding `start` will be ignored.
#'
#'See [get_params_flat()], [get_params_1comp()], [get_params_2comp()] for
#'examples.
#'
#'@param name Character: The name of the model
#'@param params Character vector: Parameter names of the model
#'@param conc_fun Name of a function to predict tissue concentrations. See
#'  Details for requirements.
#'@param auc_fun Name of the function to predict AUC (area under the
#'  concentration-time curve). See Details for requirements.
#'@param params_fun Name of the function that produces the `data.frame` of
#'  parameter info (see Details)
#'@param conc_fun_args Any additional arguments to `conc_fun` other than those
#'  listed in Details. Default NULL.
#'@param auc_fun_args Any additional arguments to `auc_fun` other those those
#'  listed in Details. Default NULL.
#'@param params_fun_args Any additional arguments to `params_fun` other than
#'  `data` (see Details). Default NULL.
#'@return An object of class `pk_model`. Effectively, a named list containing
#'  all of the arguments provided to this function.
#'@author Caroline Ring
#'@export
pk_model <- function(name,
                     params,
                     conc_fun,
                     auc_fun,
                     params_fun,
                     conc_fun_args = NULL,
                     auc_fun_args = NULL,
                     params_fun_args = NULL,
                     ...){
  #get arguments and values as a list
  argg <- c(as.list(environment()), list(...))
  this_model <- argg
  #set class
  class(this_model) <- c(class(this_model), "pk_model")
  return(this_model)
}
