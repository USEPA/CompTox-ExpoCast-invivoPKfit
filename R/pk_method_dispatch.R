#' Preprocess data generic
#' @param obj the pk object.
#' @param ... Additional arguments currently not in use.
#' @export
#' @seealso [do_preprocess.pk()] for the `do_preprocess` method for class [pk()]
do_preprocess <- function(obj, ...){
  UseMethod("do_preprocess", obj)
}

#' do_preprocess default method
#' @param obj an object
#'@param ... Additional arguments currently not in use.
#' @export
do_preprocess.default <- function(obj, ...){
  stop(paste("No 'do_preprocess' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' do_data_info generic
#' @param obj the pk object
#'@param ... Additional arguments currently not in use.
#' @export
#' @seealso [do_data_info.pk()] for the `do_data_info` method for class [pk()]
do_data_info <- function(obj, ...){
  UseMethod("do_data_info", obj)
}

#' do_data_info default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @export
do_data_info.default <- function(obj, ...){
  stop(paste("No 'do_data_info' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Prefitting
#' @param obj the pk object
#' @param ... Additional arguments currently not in use.
#' @export
do_prefit <- function(obj, ...){
  UseMethod("do_prefit", obj)
}

#' do_prefit default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#'
#' @export
do_prefit.default <- function(obj, ...){
  stop(paste("No 'do_prefit' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Fitting
#'
#' This is the S3 generic method for `do_fit`.
#' @param obj the pk object
#' @param ... Additional arguments currently not in use.
#' @export
#' @seealso [do_fit.pk()] for the `do_fit` method for class [pk()]
do_fit <- function(obj, ...){
  UseMethod("do_fit", obj)
}

#' do_fit default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @export
do_fit.default <- function(obj, ...){
  stop(paste("No 'do_fit' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Root mean squared error (RMSE)
#'
#' This is the S3 method generic for `rmse`.
#' @param obj the pk object
#' @param ... Additional arguments currently not in use.
#' @export
#' @seealso [rmse.pk()] for the `rmse` method for class [pk()]
rmse <- function(obj, ...){
  UseMethod("rmse", obj)
}

#' Root mean squared error (RMSE) default method
#' @param obj an object
#' @param ... Additional arguments.
#' @export
rmse.default <- function(obj, ...){
  stop(paste("No 'rmse' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Coefficient standard deviations
#'
#' This is the S3 method generic for `coef_sd`.
#'
#' @param obj An object
#' @param model The TK model used.
#' @param method Optimizer method used.
#' @param suppress.messages Boolean. Whether messages will be printed.
#' @param ... Additional arguments currently not in use.
#' @export
#' @seealso [coef_sd.pk()] for the `coef_sd` method for class [pk()]
coef_sd <- function(obj,
                    model,
                    method,
                    suppress.messages, ...){
  UseMethod("coef_sd", obj)
}

#' Coefficient standard deviation default
#'
#' @param obj An object
#' @param ... Additional arguments currently not in use.
#' @export
coef_sd.default <- function(obj, ...){
  stop(paste("No 'coef_sd' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Get status
#'
#' This is the S3 method generic.
#'
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @seealso [get_status.pk()] for the `get_status` method for class [pk()]
#' @export
get_status <- function(obj, ...){
  UseMethod("get_status", obj)
}

#' Default method for getting status
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @export
get_status.default <- function(obj, ...){
  stop(paste("No 'get_status' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Fold error
#'
#' This is the S3 method generic for `fold_errors`.
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @export
#' @seealso [fold_errors.pk()] for the `fold_errors` method for class [pk()]
fold_errors <- function(obj, ...){
  UseMethod("fold_errors", obj)
}

#' Fold_error default method
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @export
fold_errors.default <- function(obj, ...){
  stop(paste("No 'fold_errors' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Check required status
#'
#' This is the S3 method generic.
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [check_required_status.pk()] for the method for class [pk()]
#' @export
check_required_status <- function(obj, ...){
  UseMethod("check_required_status", obj)
}

#' Default method for checking required status
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @export
check_required_status.default <- function(obj, ...){
  stop(paste("No 'check_required_status' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Get TK stats
#'
#' This is the S3 method generic for get_tkstats(0)
#'
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @seealso [get_tkstats.pk()] for the method for class [pk()]
#' @export
get_tkstats <- function(obj, ...){
  UseMethod("get_tkstats", obj)
}

#' Default method for get_tkstats()
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @export
get_tkstats.default <- function(obj, ...){
  stop(paste("No 'get_tkstats' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Model comparison
#'
#' This is the S3 method generic for compare_models()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [compare_models.pk()] for the method for class [pk()]
#' @export
compare_models <- function(obj, ...){
  UseMethod("compare_models", obj)
}

#' Default method for compare_models()
#' @param obj an object
#' @param ... Additional arguments currently not in use.
#' @export
compare_models.default <- function(obj, ...){
  stop(paste("No 'compare_models' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}


#' get_data()
#'
#' This is the S3 method generic for get_data()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_data.pk()] for the method for class [pk()]
#' @export
get_data <- function(obj, ...){
  UseMethod("get_data", obj)
}

#' Default method for get_data()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_data.default <- function(obj, ...){
  stop(paste("No 'get_data' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_nca()
#'
#' This is the S3 method generic for get_nca()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_nca.pk()] for the method for class [pk()]
#' @export
get_nca <- function(obj, ...){
  UseMethod("get_nca", obj)
}

#' Default method for get_nca()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_nca.default <- function(obj, ...){
  stop(paste("No 'get_nca' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_data_info()
#'
#' This is the S3 method generic for get_data_info()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_data_info.pk()] for the method for class [pk()]
#' @export
get_data_info <- function(obj, ...){
  UseMethod("get_data_info", obj)
}

#' Default method for get_data_info()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_data_info.default <- function(obj, ...){
  stop(paste("No 'get_data_info' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_prefit()
#'
#' This is the S3 method generic for get_prefit()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_prefit.pk()] for the method for class [pk()]
#' @export
get_prefit <- function(obj, ...){
  UseMethod("get_prefit", obj)
}

#' Default method for get_prefit()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_prefit.default <- function(obj, ...){
  stop(paste("No 'get_prefit' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_settings_preprocess()
#'
#' This is the S3 method generic for get_settings_preprocess()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_settings_preprocess.pk()] for the method for class [pk()]
#' @export
get_settings_preprocess <- function(obj, ...){
  UseMethod("get_settings_preprocess", obj)
}

#' Default method for get_settings_preprocess()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_settings_preprocess.default <- function(obj, ...){
  stop(paste("No 'get_settings_preprocess' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_settings_data_info()
#'
#' This is the S3 method generic for get_settings_data_info()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_settings_data_info.pk()] for the method for class [pk()]
#' @export
get_settings_data_info <- function(obj, ...){
  UseMethod("get_settings_data_info", obj)
}

#' Default method for get_settings_data_info()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_settings_data_info.default <- function(obj, ...){
  stop(paste("No 'get_settings_data_info' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_stat_error_model()
#'
#' This is the S3 method generic for get_stat_error_model()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_stat_error_model.pk()] for the method for class [pk()]
#' @export
get_stat_error_model <- function(obj, ...){
  UseMethod("get_stat_error_model", obj)
}

#' Default method for get_stat_error_model()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_stat_error_model.default <- function(obj, ...){
  stop(paste("No 'get_stat_error_model' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_data_sigma_group()
#'
#' This is the S3 method generic for get_data_sigma_group()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_data_sigma_group.pk()] for the method for class [pk()]
#' @export
get_data_sigma_group <- function(obj, ...){
  UseMethod("get_data_sigma_group", obj)
}

#' Default method for get_data_sigma_group()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_data_sigma_group.default <- function(obj, ...){
  stop(paste("No 'get_data_sigma_group' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_settings_optimx()
#'
#' This is the S3 method generic for get_settings_optimx()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_settings_optimx.pk()] for the method for class [pk()]
#' @export
get_settings_optimx <- function(obj, ...){
  UseMethod("get_settings_optimx", obj)
}

#' Default method for get_settings_optimx()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_settings_optimx.default <- function(obj, ...){
  stop(paste("No 'get_settings_optimx' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_scale_conc()
#'
#' This is the S3 method generic for get_scale_conc()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_scale_conc.pk()] for the method for class [pk()]
#' @export
get_scale_conc <- function(obj, ...){
  UseMethod("get_scale_conc", obj)
}

#' Default method for get_scale_conc()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_scale_conc.default <- function(obj, ...){
  stop(paste("No 'get_scale_conc' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_scale_time()
#'
#' This is the S3 method generic for get_scale_time()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_scale_time.pk()] for the method for class [pk()]
#' @export
get_scale_time <- function(obj, ...){
  UseMethod("get_scale_time", obj)
}

#' Default method for get_scale_time()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_scale_time.default <- function(obj, ...){
  stop(paste("No 'get_scale_time' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_error_group()
#'
#' This is the S3 method generic for get_error_group()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_error_group.pk()] for the method for class [pk()]
#' @export
get_error_group <- function(obj, ...){
  UseMethod("get_error_group", obj)
}

#' Default method for get_error_group()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_error_group.default <- function(obj, ...){
  stop(paste("No 'get_error_group' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_stat_model()
#'
#' This is the S3 method generic for get_stat_model()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_stat_model.pk()] for the method for class [pk()]
#' @export
get_stat_model <- function(obj, ...){
  UseMethod("get_stat_model", obj)
}

#' Default method for get_stat_model()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_stat_model.default <- function(obj, ...){
  stop(paste("No 'get_stat_model' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_data_original()
#'
#' This is the S3 method generic for get_data_original()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_data_original.pk()] for the method for class [pk()]
#' @export
get_data_original <- function(obj, ...){
  UseMethod("get_data_original", obj)
}

#' Default method for get_data_original()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_data_original.default <- function(obj, ...){
  stop(paste("No 'get_data_original' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_mapping()
#'
#' This is the S3 method generic for get_mapping()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_mapping.pk()] for the method for class [pk()]
#' @export
get_mapping <- function(obj, ...){
  UseMethod("get_mapping", obj)
}

#' Default method for get_mapping()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_mapping.default <- function(obj, ...){
  stop(paste("No 'get_mapping' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' rsq()
#'
#' This is the S3 method generic for rsq()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [rsq.pk()] for the method for class [pk()]
#' @export
rsq <- function(obj, ...){
  UseMethod("rsq", obj)
}

#' Default method for rsq()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
rsq.default <- function(obj, ...){
  stop(paste("No 'rsq' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_winning_model()
#'
#' This is the S3 method generic for get_winning_model()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [get_winning_model.pk()] for the method for class [pk()]
#' @export
get_winning_model <- function(obj, ...){
  UseMethod("get_winning_model", obj)
}

#' Default method for get_winning_model()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_winning_model.default <- function(obj, ...){
  stop(paste("No 'get_winning_model' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' nca()
#'
#' This is the S3 method generic for nca()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [nca.pk()] for the method for class [pk()]
#' @export
nca <- function(obj, ...){
  UseMethod("nca", obj)
}

#' Default method for nca()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
nca.default <- function(obj, ...){
  stop(paste("No 'nca' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' data_summary()
#'
#' This is the S3 method generic for data_summary()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
data_summary <- function(obj, ...){
  UseMethod("data_summary", obj)
}

#' Default method for data_summary()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
data_summary.default <- function(obj, ...){
  stop(paste("No 'data_summary' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}


#' get_data_summary()
#'
#' This is the S3 method generic for get_data_summary()
#'
#' `get_data_summary()` is an alias for `data_summary()`
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
get_data_summary <- function(obj, ...){
  UseMethod("data_summary", obj)
}

#' Default method for get_data_summary()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_data_summary.default <- function(obj, ...){
  stop(paste("No 'get_data_summary' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' eval_tkstats()
#'
#' This is the S3 method generic for eval_tkstats()
#'
#' `eval_tkstats()` is an alias for `data_summary()`
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
eval_tkstats <- function(obj, ...){
  UseMethod("eval_tkstats", obj)
}

#' Default method for eval_tkstats()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
eval_tkstats.default <- function(obj, ...){
  stop(paste("No 'eval_tkstats' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_fit()
#'
#' This is the S3 method generic for get_fit()
#'
#' `get_fit()` is an alias for `data_summary()`
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
get_fit <- function(obj, ...){
  UseMethod("get_fit", obj)
}

#' Default method for get_fit()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_fit.default <- function(obj, ...){
  stop(paste("No 'get_fit' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_data_group()
#'
#' This is the S3 method generic for get_data_group()
#'
#' `get_data_group()` is an alias for `data_summary()`
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
get_data_group <- function(obj, ...){
  UseMethod("get_data_group", obj)
}

#' Default method for get_data_group()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
get_data_group.default <- function(obj, ...){
  stop(paste("No 'get_data_group' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' twofold_test()
#'
#' This is the S3 method generic for twofold_test()
#'
#' @param obj An object.
#' @param ... Additional arguments currently not in use.
#' @seealso [twofold_test.pk()] for the method for class [pk()]
#' @export
twofold_test <- function(obj, ...){
  UseMethod("twofold_test", obj)
}

#' Default method for twofold_test()
#'
#'@param obj An object
#'@param ... Additional arguments currently not in use.
#' @export
twofold_test.default <- function(obj, ...){
  stop(paste("No 'twofold_test' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}
