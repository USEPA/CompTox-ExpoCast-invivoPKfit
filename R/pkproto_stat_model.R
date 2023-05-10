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

#'Error model
#'
#'Define an error model.
#'
#'`stat_error_model` defines groupings for a fixed-effects error model. For each
#'model in `stat_model`, a single set of model parameters will be fit to `data`.
#'In order to do the fitting, the residual errors (observed concentrations -
#'model-predicted concentrations) are assumed to obey a zero-mean normal
#'distribution. However, in this package, the residuals are not all required to
#'obey the *same* zero-mean normal distribution. Different groups of residuals
#'may obey zero-mean normal distributions with different variances.
#'`stat_error_model` defines these groups as unique combinations of the
#'variables given in argument `error_group`. For example, the default value
#'`vars(Chemical, Species, Reference, Media)` means that for each group of
#'observations in `data` with a unique combination of `Chemical`, `Species`,
#'`Reference`, and `Media`, there is a separate residual error variance. For
#'example, if there happened to be three such unique combinations, there would
#'be three error variances.
#'
#'If you want all residuals to obey the same zero-mean normal distribution
#'(i.e., for there to be only one residual error variance), then you should
#'provide an `error_group` that puts all the data in the same group. For
#'example, since all data in `data` should already be for a single `Chemical`
#'and `Species`, you could provide `error_group = vars(Chemical, Species)` to
#'put all the data in the same group.
#'
#'Note that, since all data in `data` should already be for a single `Chemical`
#'and `Species`, you could leave out `Chemical` and `Species` from `error_group`
#'and still get the same result. However, we recommend explicitly including
#'`Chemical` and `Species`. Tncluding them will make your code more explicit and
#'transparent, and it does no harm. In addition, [invivopkfit] may be extended
#'in the future to allow input of data with multiple chemicals or species;
#'explicitly including `Chemical` and `Species` in your `error_group` will
#'future-proof your code in that sense.
#'
#'The error variance(s) are hyperparameters that will be estimated from the data
#'along with the model parameters. That means there needs to be enough data to
#'fit the model parameters plus the error variances. For example, if you are
#'fitting a 1-compartment model to oral and IV data measured in plasma, and
#'using an error model with three separate error-variance groups (e.g. three
#'different References), then you are trying to fit 4 model parameters (kelim,
#'Vdist, Fgutabs, kgutabs) plus 3 error variances, for a total of 7 parameters.
#'That means you need to have at least 8 data points. (When you call [prefit()],
#'this checking is done automatically. But it is useful to be aware of this, in
#'case you are trying to figure out why your fit was aborted due to insufficient data availability.)
#'
#'@param error_group Defined using [dplyr::vars()]: A set of harmonized
#'  variables whose unique combinations define a group with its own error
#'  variance. These variables refer to the `data` element of the `pk` object to
#'  which `stat_error_model()` is added. The variables should not be quoted.
#'  Default is `vars(Chemical, Species, Reference, Media)`.
#'@param ... Additional arguments. Not currently used.
#'@return An object of class `pk_stat_error_model`: A named list of all the
#'  arguments to `stat_error_model`.
#'@author Caroline Ring
#'@export
stat_error_model <- function(error_group = vars(Chemical, Species, Reference, Media),
                             ...){
  #get arguments and values as a list
  argg <- c(as.list(environment()), list(...))
  this_error_model <- argg
  class(this_error_model) <- c(class(this_error_model), "pkproto", "pk_stat_error_model")
  return(this_error_model)
}

# stat_nca_model <- function(nca_group = vars(Chemical, Species, Reference, Route, Media),
#                            )
