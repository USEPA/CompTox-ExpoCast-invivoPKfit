#' Add PK model(s) to be fitted
#'
#' @param model The name(s) of models to be fitted. These should be the names of objects of class `pk_model`. Built-in options are `"flat"`,
#'   `"1comp"`, and `"2comp"`. You may add your own model by using [pk_model()].

stat_model.pk <- function(model = c("flat", "1comp", "2comp"),
                      ...){
  this_stat_model <- list()
  #pull pk_model objects by each name
for(this_model in model){
  this_stat_model[[this_model]] <- list()
  #check whether an object exists by this name
  if(exists(this_model)){
    this_model_obj <- eval(parse(text = this_model))
    #Check whether this is an object of class `pk_model`
    if(inherits(this_model_obj, "pk_model")){
      #if so, add it to the list
      this_stat_model[[this_model]] <- this_model_obj
    }
  }

  return(this_stat_model)
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
#'# `get_params_fun` requirements
#'
#'`get_params_fun` should be a function whose first argument is a `data.frame`,
#'which will be the pre-processed data using `invivopkfit` harmonized variable
#'names. It should then return a `data.frame` with the following variables:
#'
#' - `param_name`: Character vector, listing parameter names for the model
#' - `param_units`: Character vector, listing units of each model parameter
#' - `optimize_param`: Logical (TRUE/FALSE), whether each parameter is to be estimated given the available data
#' - `use_param`: Logical (TRUE/FALSE), whether each parameter is to be used in the model even if it is not estimated (i.e., if a parameter value is to be held constant while the others are estimated, then `optimize_param` should be FALSE but `use_param` should be TRUE)
#' -`lower_bound`: Numerical. Lower bounds for each parameter. May be `-Inf` if no lower bound.  If `opt_param` or `use_param` is FALSE, thenthe corresponding `lower_bound` will be ignored (because the parameter is not being estimated from the data).
#' - `upper_bound`: Numerical. Upper bounds for each parameter. May be `Inf` if no upper bound. If `opt_param` or `use_param` is FALSE, then the corresponding `upper_bound` will be ignored (because the parameter is not being estimated from the data).
#' - `start`: Numerical. Starting values for estimating each parameter. If `optimize_param` is FALSE and `use_param` is TRUE, then the parameter will be held constant at the corresponding value in `start`. If `use_param` is FALSE, then the corresponding `start` will be ignored.
#'
#' See [get_params_flat()], [get_params_1comp()], [get_params_2comp()] for examples.
#'
#'
#'
#'@param name Character: The name of the model
#'@param params Character vector: Parameter names of the model
#'@param conc_fun Name of a function to predict tissue concentrations. See
#'  Details for requirements.
#'@param auc_fun Name of the function to predict AUC (area under the
#'  concentration-time curve). See Details for requirements.
#'@param params_fun Name of the function that produces the `data.frame` of parameter info (see Details)
#' @param conc_fun_args Any additional arguments to `conc_fun` other than those listed in Details. Default NULL.
#' @param auc_fun_args Any additional arguments to `auc_fun` other those those listed in Details. Default NULL.
#' @param params_fun_args Any additional arguments to `params_fun` other than `data` (see Details). Default NULL.
#' @return An object of class `pk_model`. Effectively, a named list containing all of the arguments provided to this function.
#' @author Caroline Ring
#' @export
pk_model <- function(name,
                     params,
                     conc_fun,
                     auc_fun,
                     params_fun,
                     conc_fun_args = NULL,
                     auc_fun_args = NULL,
                     params_fun_args = NULL){
  #get arguments and values as a list
  argg <- c(as.list(environment()), list(...))
  this_model <- argg
  #set class
  class(this_model) <- c(class(this_model), "pk_model")
  return(this_model)
}

#'Error model
#'
#'
#'@param error_model One of "FE" (the default) or "pooled". "FE" stands for
#'  "fixed-effects" and means that there is only one set of model parameters,
#'  but each group defined in `error_group` will have its own error variance
#'  around the estimated concentration-time curve. "pooled" means that all data
#'  is assumed to share the same error variance. If `error_model ='pooled'` then
#'  `error_group` will be ignored.
#'@param error_group Defined using [ggplot2::vars()]: A set of harmonized
#'  variables whose unique combinations define a group with its own error
#'  variance. If `error_group`is chosen such that all data are in the same group
#'  (i.e. there is only one group), then the effective error model will be the
#'  same regardless of `error_model`.
#'@param ... Additional arguments. Not currently used.
#'@return An object of class `pk_stat_error_model`: A named list of all the
#'  arguments to `stat_error_model`.
#'@author Caroline Ring
#'@export
stat_error_model <- function(error_model ="FE",
                             error_group = vars(Chemical, Species, Reference, Media),
                             ...){
  #get arguments and values as a list
  argg <- c(as.list(environment()), list(...))
  this_error_model <- argg
  class(this_error_model) <- c(class(this_model), "pk_stat_error_model")
  return(this_error_model)
}
