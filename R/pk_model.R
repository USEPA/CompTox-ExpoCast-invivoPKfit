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
#'# `tkstats_fun` requirements
#'
#'`tkstats_fun` should be a function which accepts a vector of model parameter
#'values and calculates derived summary toxicokinetic statistics (e.g. total
#'clearance, halflife, AUC, volume of distribution at steadystate).
#'
#'The function must take the following named arguments:
#'
#' - `pars`: A named numeric vector of model parameter values.
#' - `route`: A character scalar naming a route (e.g. "oral" or "iv")
#' - `medium`: A character scalar naming a tissue medium of analysis (e.g. "blood" or "plasma")
#' - `dose`: A numeric scalar giving a dose level for which to calculate TK statistics
#' - `time_unit`: A character scalar giving the units of time
#' - `conc_unit`: A character scalar giving the units of concentration
#' - `vol_unit`: A character scalar giving the units of volume
#'
#'and return a `data.frame` of derived toxicokinetic statistics, which should
#'have the following variables:
#'
#' - `param_name`: A character vector giving the names of each derived TK statistic
#' - `param_value`: A character vector giving the values of each derived TK statistic
#' - `param_units`: A character vector giving the units of each derived TK statistic
#'
#'It is recommended (although not required) that the function return the
#'following statistics, using these names in the `param_name` variable:
#'
#' - `CLtot`: Total clearance rate (units of volume/time)
#' - `CLtot/Fgutabs`: Total clearance rate scaled by bioavailability (if oral bioavailability is available) (units of volume/time)
#' - `Css`: The steady-state concentration after a long-term daily dose of `dose` (units of concentration)
#' - `halflife`: The terminal half-life (units of time)
#' - `tmax`: The time of peak concentration (units of time)
#' - `Cmax`: The peak concentration (units of time)
#' - `AUC_infinity`: The area under the concentration-time curve, calculated at infinite time (units of concentration * time)
#' - `Vdist_ss`: The volume of distribution at steady-state (units of volume)
#' - `Vdist_ss/Fgutabs`: The volume of distribution at steady-state scaled by bioavailability (if oral bioavailability is available) (units of volume)
#'
#'The recommendation to return these statistics, using these names, is intended
#'to make it easier to compare TK statistics across models, and to compare TK
#'statistics to the results of non-compartmental analysis. If these names are
#'not used, then some outputs of [summary.pk()] will not be very useful. The
#'automated comparison of TK stats from the winning model to the results of
#'non-compartmental analysis relies on these names being present in the output
#'of `tkstats_fun` to match the names of the statistics output from NCA; it
#'shouldn't crash without them, but the results won't be very useful. And TK
#'stats compiled across models will not be easy to compare if the models use
#'different names for the statistics.
#'
#'@param name Character: The name of the model.
#'@param params Character vector: Parameter names of the model.
#'@param conc_fun Character: Name of the function to predict tissue concentrations using this model. See
#'  Details for requirements.
#'@param auc_fun Character: Name of the function to predict AUC (area under the
#'  concentration-time curve) using this model. See Details for requirements.
#'@param params_fun Character: Name of the function that produces the `data.frame` of
#'  parameter info for this model (see Details)
#'@param tkstats_fun Character: Name of the function that produces a`data.frame` of derived
#'  TK statistics for this model (see Details)
#'@param conc_fun_args A named list: any additional arguments to `conc_fun` other than those
#'  listed in Details. Default `NULL`.
#'@param auc_fun_args A named list: any additional arguments to `auc_fun` other those those
#'  listed in Details. Default `NULL`.
#'@param params_fun_args A named list: any additional arguments to `params_fun` other than
#'  `data` (see Details). Default `NULL`.
#'@param tkstats_fun_args A named list: any additional arguments to `tkstats_fun` other than
#'  `data`, `medium`, `route` (see Details). Default NULL.
#'@param ... Additional arguments (not currently implemented).
#'@return An object of class `pk_model`. Effectively, a named list containing
#'  all of the arguments provided to this function.
#'@author Caroline Ring
#'@export
pk_model <- function(name,
                     params,
                     conc_fun,
                     auc_fun,
                     params_fun,
                     tkstats_fun,
                     conc_fun_args = NULL,
                     auc_fun_args = NULL,
                     params_fun_args = NULL,
                     tkstats_fun_args = NULL,
                     ...){
  #get arguments and values as a list
  argg <- c(as.list(environment()), list(...))
  this_model <- argg
  #set class
  class(this_model) <- c(class(this_model), "pk_model")
  return(this_model)
}
