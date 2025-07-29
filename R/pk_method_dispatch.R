#' Preprocess data generic
#' @param obj the pk object.
#' @param ... Additional arguments currently not in use.
#' @return The same `pk` object, with added elements `data` (containing the
#'   cleaned, gap-filled data) and `data_info` (containing summary information
#'   about the data, e.g. number of observations by route, media,
#'   detect/nondetect; empirical tmax, time of peak concentration for oral data;
#'   number of observations before and after empirical tmax)
#' @export
#' @seealso [do_preprocess.pk()] for the `do_preprocess` method for class [pk()]
do_preprocess <- function(obj, ...) {
  UseMethod("do_preprocess", obj)
}

#' do_preprocess default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
do_preprocess.default <- function(obj, ...) {
  stop("No 'do_preprocess' method exists for object of class",
       toString(class(obj))
  )
}

#' do_data_info generic
#' @param obj the pk object
#' @param ... Additional arguments currently not in use.
#' @return Object of class [pk()] with an added `$data_info` list containing
#' non-compartmental analysis results.
#' @export
#' @seealso [do_data_info.pk()] for the `do_data_info` method for class [pk()]
do_data_info <- function(obj, ...) {
  UseMethod("do_data_info", obj)
}

#' do_data_info default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
do_data_info.default <- function(obj, ...) {
  stop("No 'do_data_info' method exists for object of class",
       toString(class(obj))
  )
}

#' Prefitting
#' @param obj the pk object
#' @param ... Additional arguments currently not in use.
#' @return The same `pk` object, but with a new element `prefit`, containing the
#'   results of pre-fit calculations and checks for each model and for the error
#'   model.
#' @export
#' @seealso [do_prefit.pk()] for the `do_prefit` method for class [pk()]
do_prefit <- function(obj, ...) {
  UseMethod("do_prefit", obj)
}

#' do_prefit default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#'
#' @export
do_prefit.default <- function(obj, ...) {
  stop("No 'do_prefit' method exists for object of class",
       toString(class(obj))
  )
}

#' Fitting
#'
#' This is the S3 generic method for `do_fit`.
#' @param obj the pk object
#' @param ... Additional arguments currently not in use.
#' @return The same [pk] object, with element `fit` containing the fitted
#'   results for each model in `stat_model`.
#' @export
#' @seealso [do_fit.pk()] for the `do_fit` method for class [pk()]
do_fit <- function(obj, ...) {
  UseMethod("do_fit", obj)
}

#' do_fit default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
do_fit.default <- function(obj, ...) {
  stop("No 'do_fit' method exists for object of class",
       toString(class(obj))
  )
}

#' Root mean squared error (RMSE)
#'
#' This is the S3 method generic for `rmse`.
#' @param obj the pk object
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame` with calculated RMSE as the final column. There is one row per
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data, each [optimx::optimx()] methods (specified in
#'   [settings_optimx()]), `rmse_group` specified.
#' @export
#' @seealso [rmse.pk()] for the `rmse` method for class [pk()]
rmse <- function(obj, ...) {
  UseMethod("rmse", obj)
}

#' Root mean squared error (RMSE) default method
#' @param obj an object
#' @param ... Additional arguments.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
rmse.default <- function(obj, ...) {
  stop("No 'rmse' method exists for object of class",
       toString(class(obj))
  )
}


#' Coefficient standard deviations
#'
#' This is the S3 method generic for `coef_sd`.
#'
#' @param obj A pk object.
#' @param model The TK model used.
#' @param method Optimizer method used.
#' @param suppress.messages Boolean. Whether messages will be printed.
#' @param ... Additional arguments currently not in use.
#' @return A dataframe with one row for each `data_group`, `model` and `method`.
#'   The remaining columns include the parameters & hyperparameters as returned
#'   by [coef.pk()], as well as their calculated standard deviations.
#' @export
#' @seealso [coef_sd.pk()] for the `coef_sd` method for class [pk()]
coef_sd <- function(obj,
                    model,
                    method,
                    suppress.messages, ...) {
  UseMethod("coef_sd", obj)
}

#' Coefficient standard deviation default
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
coef_sd.default <- function(obj, ...) {
  stop("No 'coef_sd' method exists for object of class",
       toString(class(obj))
  )
}

#' Get status
#'
#' This is the S3 method generic.
#'
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return The status of the `pk` object as an integer.
#' @seealso [get_status.pk()] for the `get_status` method for class [pk()]
#' @export
get_status <- function(obj, ...) {
  UseMethod("get_status", obj)
}

#' Default method for getting status
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_status.default <- function(obj, ...) {
  stop("No 'get_status' method exists for object of class",
       toString(class(obj))
  )
}

#' Fold error
#'
#' This is the S3 method generic for `fold_error`.
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return  A data.frame with one row for each `data_group`, `model` and
#'   `method`. A column contains the fold errors (observed/predicted) of the
#'   model fitted by the corresponding method. These residuals are
#'   concentrations in the same units as `obj$data$Conc.Units`; any
#'   concentration transformations (in `obj$scale$conc`) are *not* applied.
#' @export
#' @seealso [fold_error.pk()] for the `fold_error` method for class [pk()]
fold_error <- function(obj, ...) {
  UseMethod("fold_error", obj)
}

#' fold_error default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
fold_error.default <- function(obj, ...) {
  stop("No 'fold_error' method exists for object of class",
       toString(class(obj))
  )
}


#' Check required status
#'
#' This is the S3 method generic.
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return If the [pk()] object has the required status or greater, returns
#'   TRUE. If the [pk()] object has less than the required status, returns
#'   FALSE. Returned value has an attribute `msg`, containing an informative
#'   message as a string.
#' @seealso [check_required_status.pk()] for the method for class [pk()]
#' @export
check_required_status <- function(obj, ...) {
  UseMethod("check_required_status", obj)
}

#' Default method for checking required status
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
check_required_status.default <- function(obj, ...) {
  stop("No 'check_required_status' method exists for object of class",
       toString(class(obj))
  )
}

#' Get TK stats
#'
#' This is the S3 method generic for get_tkstats(0)
#'
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return  A data.frame with one row for each `data_group`, `model` and `method`
#'   with the variables in the `data.frame` returned by the `tkstats_fun` for
#'   its corresponding model.
#'   (For the built-in models `model_flat`, `model_1comp`, and `model_2comp`, these
#'   variables are `param_name` and `param_value`.)
#' @seealso [get_tkstats.pk()] for the method for class [pk()]
#' @export
get_tkstats <- function(obj, ...) {
  UseMethod("get_tkstats", obj)
}

#' Default method for get_tkstats()
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_tkstats.default <- function(obj, ...) {
  stop("No 'get_tkstats' method exists for object of class",
       toString(class(obj))
  )
}

#' Model comparison
#'
#' This is the S3 method generic for compare_models()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame` with variables
#' - `model`: The name of each model
#' - `method`: The name of each method
#' - A variable named for `criterion` (e.g. if `criterion = "AIC"` then the result will have a
#'   variable named `AIC`): The criterion value for each model/method
#' @seealso [compare_models.pk()] for the method for class [pk()]
#' @export
compare_models <- function(obj, ...) {
  UseMethod("compare_models", obj)
}

#' Default method for compare_models()
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
compare_models.default <- function(obj, ...) {
  stop("No 'compare_models' method exists for object of class",
       toString(class(obj))
  )
}


#' get_data()
#'
#' This is the S3 method generic for get_data()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame`: the `data` element of `obj`
#' @seealso [get_data.pk()] for the method for class [pk()]
#' @export
get_data <- function(obj, ...) {
  UseMethod("get_data", obj)
}

#' Default method for get_data()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_data.default <- function(obj, ...) {
  stop("No 'get_data' method exists for object of class",
       toString(class(obj))
  )
}

#' get_nca()
#'
#' This is the S3 method generic for get_nca()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame`: the `data` element of `obj`
#' @seealso [get_nca.pk()] for the method for class [pk()]
#' @export
get_nca <- function(obj, ...) {
  UseMethod("get_nca", obj)
}

#' Default method for get_nca()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_nca.default <- function(obj, ...) {
  stop("No 'get_nca' method exists for object of class",
       toString(class(obj))
  )
}

#' get_data_info()
#'
#' This is the S3 method generic for get_data_info()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `list` of `tibble`s: the `data_info` element of `obj`
#' @seealso [get_data_info.pk()] for the method for class [pk()]
#' @export
get_data_info <- function(obj, ...) {
  UseMethod("get_data_info", obj)
}

#' Default method for get_data_info()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_data_info.default <- function(obj, ...) {
  stop("No 'get_data_info' method exists for object of class",
       toString(class(obj))
  )
}

#' get_prefit()
#'
#' This is the S3 method generic for get_prefit()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_prefit.pk()] for the method for class [pk()]
#' @return A list of `data.frame`s in the object's `prefit` element.
#' @export
get_prefit <- function(obj, ...) {
  UseMethod("get_prefit", obj)
}

#' Default method for get_prefit()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_prefit.default <- function(obj, ...) {
  stop("No 'get_prefit' method exists for object of class",
       toString(class(obj))
  )
}

#' get_settings_preprocess()
#'
#' This is the S3 method generic for get_settings_preprocess()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A named list of the preprocessing settings
#' @seealso [get_settings_preprocess.pk()] for the method for class [pk()]
#' @export
get_settings_preprocess <- function(obj, ...) {
  UseMethod("get_settings_preprocess", obj)
}

#' Default method for get_settings_preprocess()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_settings_preprocess.default <- function(obj, ...) {
  stop("No 'get_settings_preprocess' method exists for object of class",
       toString(class(obj))
  )
}

#' get_nca_group()
#'
#' This is the S3 method generic for get_nca_group()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A named list of the data_info settings
#' @seealso [get_nca_group.pk()] for the method for class [pk()]
#' @export
get_nca_group <- function(obj, ...) {
  UseMethod("get_nca_group", obj)
}

#' Default method for get_nca_group()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_nca_group.default <- function(obj, ...) {
  stop("No 'get_nca_group' method exists for object of class",
       toString(class(obj))
  )
}

#' get_data_sigma_group()
#'
#' This is the S3 method generic for get_data_sigma_group()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `factor` vector giving the error SD group ID for each observation,
#'   as the interaction of the factors specified in
#'   `obj$pk_groups$error_group`.
#' @seealso [get_data_sigma_group.pk()] for the method for class [pk()]
#' @export
get_data_sigma_group <- function(obj, ...) {
  UseMethod("get_data_sigma_group", obj)
}

#' Default method for get_data_sigma_group()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_data_sigma_group.default <- function(obj, ...) {
  stop("No 'get_data_sigma_group' method exists for object of class",
       toString(class(obj))
  )
}

#' get_settings_optimx()
#'
#' This is the S3 method generic for get_settings_optimx()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A named list of the optimx settings
#' @seealso [get_settings_optimx.pk()] for the method for class [pk()]
#' @export
get_settings_optimx <- function(obj, ...) {
  UseMethod("get_settings_optimx", obj)
}

#' Default method for get_settings_optimx()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_settings_optimx.default <- function(obj, ...) {
  stop("No 'get_settings_optimx' method exists for object of class",
       toString(class(obj))
  )
}

#' get_scale_conc()
#'
#' This is the S3 method generic for get_scale_conc()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `list`: `obj$scales$conc`
#' @seealso [get_scale_conc.pk()] for the method for class [pk()]
#' @export
get_scale_conc <- function(obj, ...) {
  UseMethod("get_scale_conc", obj)
}

#' Default method for get_scale_conc()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_scale_conc.default <- function(obj, ...) {
  stop("No 'get_scale_conc' method exists for object of class",
       toString(class(obj))
  )
}

#' get_scale_time()
#'
#' This is the S3 method generic for get_scale_time()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `list`: `obj$scales$time`
#' @seealso [get_scale_time.pk()] for the method for class [pk()]
#' @export
get_scale_time <- function(obj, ...) {
  UseMethod("get_scale_time", obj)
}

#' Default method for get_scale_time()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_scale_time.default <- function(obj, ...) {
  stop("No 'get_scale_time' method exists for object of class",
       toString(class(obj))
  )
}

#' get_error_group()
#'
#' This is the S3 method generic for get_error_group()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return The stat_error_model error grouping
#' @seealso [get_error_group.pk()] for the method for class [pk()]
#' @export
get_error_group <- function(obj, ...) {
  UseMethod("get_error_group", obj)
}

#' Default method for get_error_group()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_error_group.default <- function(obj, ...) {
  stop("No 'get_error_group' method exists for object of class",
       toString(class(obj))
  )
}

#' get_stat_model()
#'
#' This is the S3 method generic for get_stat_model()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `list` -- the `stat_model` element of `obj`
#' @seealso [get_stat_model.pk()] for the method for class [pk()]
#' @export
get_stat_model <- function(obj, ...) {
  UseMethod("get_stat_model", obj)
}

#' Default method for get_stat_model()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_stat_model.default <- function(obj, ...) {
  stop("No 'get_stat_model' method exists for object of class",
       toString(class(obj))
  )
}

#' get_data_original()
#'
#' This is the S3 method generic for get_data_original()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame` -- the `data_original` element of `obj`
#' @seealso [get_data_original.pk()] for the method for class [pk()]
#' @export
get_data_original <- function(obj, ...) {
  UseMethod("get_data_original", obj)
}

#' Default method for get_data_original()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_data_original.default <- function(obj, ...) {
  stop("No 'get_data_original' method exists for object of class",
       toString(class(obj))
  )
}

#' get_mapping()
#'
#' This is the S3 method generic for get_mapping()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A list of `quosure`s -- the `mapping` element of `obj`
#' @seealso [get_mapping.pk()] for the method for class [pk()]
#' @export
get_mapping <- function(obj, ...) {
  UseMethod("get_mapping", obj)
}

#' Default method for get_mapping()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_mapping.default <- function(obj, ...) {
  stop("No 'get_mapping' method exists for object of class",
       toString(class(obj))
  )
}

#' rsq()
#'
#' This is the S3 method generic for rsq()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return  A dataframe with one row for each `data_group`, `model` and `method`.
#'   The final column contains the R-squared of the model fitted by the corresponding
#'   method, using the data in `newdata`.
#' @seealso [rsq.pk()] for the method for class [pk()]
#' @export
rsq <- function(obj, ...) {
  UseMethod("rsq", obj)
}

#' Default method for rsq()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
rsq.default <- function(obj, ...) {
  stop("No 'rsq' method exists for object of class",
       toString(class(obj))
  )
}

#' get_winning_model()
#'
#' This is the S3 method generic for get_winning_model()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A data.frame with one row for each `data_group`, `model` and `method` and
#'  The return value has attribute `criterion` giving the name of the criterion function used to compare
#'  models.
#' @seealso [get_winning_model.pk()] for the method for class [pk()]
#' @export
get_winning_model <- function(obj, ...) {
  UseMethod("get_winning_model", obj)
}

#' Default method for get_winning_model()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_winning_model.default <- function(obj, ...) {
  stop("No 'get_winning_model' method exists for object of class",
       toString(class(obj))
  )
}

#' nca()
#'
#' This is the S3 method generic for nca()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame` with variables including all the grouping variables in
#'   `nca_group`, `nca_group_id`; `design` (the auto-detected study design for
#'   this group); `param_name` (the name of the NCA parameter); `param_value`
#'   (the NCA parameter value); `param_sd_z` (standard deviation of the
#'   estimated NCA parameter value, if available); `param_units` (the units of
#'   the NCA parameter, derived from the units of the data).
#' @seealso [nca.pk()] for the method for class [pk()]
#' @export
nca <- function(obj, ...) {
  UseMethod("nca", obj)
}

#' Default method for nca()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
nca.default <- function(obj, ...) {
  stop("No 'nca' method exists for object of class",
       toString(class(obj))
  )
}

#' data_summary()
#'
#' This is the S3 method generic for data_summary()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame` with variables including all the grouping variables in
#'   `summary_group`, `group_id`; `param_name` (the name of the summary
#'   statistic; see Details); `param_value` (the summary statistic value);  `param_units`
#'   (the units of the summary statistic, derived from the units of the data).
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
data_summary <- function(obj, ...) {
  UseMethod("data_summary", obj)
}

#' Default method for data_summary()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
data_summary.default <- function(obj, ...) {
  stop("No 'data_summary' method exists for object of class",
       toString(class(obj))
  )
}


#' get_data_summary()
#'
#' This is the S3 method generic for get_data_summary()
#'
#' `get_data_summary()` is an alias for `data_summary()`
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame` with variables including all the grouping variables in
#'   `summary_group`, `group_id`; `param_name` (the name of the summary
#'   statistic; see Details); `param_value` (the summary statistic value);  `param_units`
#'   (the units of the summary statistic, derived from the units of the data).
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
get_data_summary <- function(obj, ...) {
  UseMethod("data_summary", obj)
}

#' Default method for get_data_summary()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_data_summary.default <- function(obj, ...) {
  stop("No 'get_data_summary' method exists for object of class",
       toString(class(obj))
  )
}

#' eval_tkstats()
#'
#' This is the S3 method generic for eval_tkstats()
#'
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A `data.frame` with one  row for each "winning" model in
#'   `model` from [get_winning_model()]. The `data.frame` will have the variables
#'   returned by the `tkstats_fun` for its corresponding model. (For the
#'   built-in models `model_flat`, `model_1comp`, and `model_2comp`, these
#'   variables are `param_name` and `param_value`.) Additionally, there will be
#'   a variable `method` denoting the [optimx::optimx()] method used to optimize
#'   the set of model parameters used to derive each set of TK statistics.
#' @seealso [eval_tkstats.pk()] for the method for class [pk()]
#' @export
eval_tkstats <- function(obj, ...) {
  UseMethod("eval_tkstats", obj)
}

#' Default method for eval_tkstats()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
eval_tkstats.default <- function(obj, ...) {
  stop("No 'eval_tkstats' method exists for object of class",
       toString(class(obj))
  )
}

#' get_fit()
#'
#' This is the S3 method generic for get_fit()
#'
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A named list of objects of class `optimx`, named for the models in
#'   `model`. As described in [optimx::optimx()]  If only one model is
#'   specified, the return value will still be a list, but with only one
#'   element.
#' @seealso [get_fit.pk()] for the method for class [pk()]
#' @export
get_fit <- function(obj, ...) {
  UseMethod("get_fit", obj)
}

#' Default method for get_fit()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_fit.default <- function(obj, ...) {
  stop("No 'get_fit' method exists for object of class",
       toString(class(obj))
  )
}

#' get_data_group()
#'
#' This is the S3 method generic for get_data_group()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return An object of class `call` giving the data grouping as a `dplyr::vars()` specification
#' @seealso [get_data_group.pk()] for the method for class [pk()]
#' @export
get_data_group <- function(obj, ...) {
  UseMethod("get_data_group", obj)
}

#' Default method for get_data_group()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_data_group.default <- function(obj, ...) {
  stop("No 'get_data_group' method exists for object of class",
       toString(class(obj))
  )
}

#' twofold_test()
#'
#' This is the S3 method generic for twofold_test()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A list of data frames.
#' @seealso [twofold_test.pk()] for the method for class [pk()]
#' @export
twofold_test <- function(obj, ...) {
  UseMethod("twofold_test", obj)
}

#' Default method for twofold_test()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
twofold_test.default <- function(obj, ...) {
  stop("No 'twofold_test' method exists for object of class",
       toString(class(obj))
  )
}

#' AFE()
#'
#' This is the S3 method generic for AFE()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return  A dataframe with one row for each `data_group`, `model` and `method`.
#'   The final column contains the AFE of the model fitted by the corresponding
#'   method, using the data in `newdata`.
#' @seealso [AFE.pk()] for the method for class [pk()]
#' @export
AFE <- function(obj, ...) {
  UseMethod("AFE", obj)
}

#' Default method for AFE()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
AFE.default <- function(obj, ...) {
  stop("No 'AFE' method exists for object of class",
       toString(class(obj))
  )
}

#' AAFE()
#'
#' This is the S3 method generic for AAFE()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return  A dataframe with one row for each `data_group`, `model` and `method`.
#'   The final column contains the AAFE of the model fitted by the corresponding
#'   method, using the data in `newdata`.
#' @seealso [AAFE.pk()] for the method for class [pk()]
#' @export
AAFE <- function(obj, ...) {
  UseMethod("AAFE", obj)
}

#' Default method for AAFE()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
AAFE.default <- function(obj, ...) {
  stop("No 'AAFE' method exists for object of class",
       toString(class(obj))
  )
}

#' get_hessian()
#'
#' This is the S3 method generic for get_hessian()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @return A dataframe with one row for each `data_group`, `model` and `method`.
#'  The remaining column is a `list` column containing the Hessian for each row.
#' @seealso [hessian.pk()] for the method for class [pk()]
#' @export
get_hessian <- function(obj, ...) {
  UseMethod("get_hessian", obj)
}

#' Default method for get_hessian()
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @return An error, when a non-pk object is used for the first argument.
#' @export
get_hessian.default <- function(obj, ...) {
  stop("No 'get_hessian' method exists for object of class",
       toString(class(obj))
  )
}



