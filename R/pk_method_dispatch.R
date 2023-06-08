#' Add a `pkproto` object to a `pk` object
#'
#' This is the S3 generic method.
#'
#' @param pkproto_obj The `pkproto` object to be added
#' @param pk_obj The `pk` object to which the `pkproto` object is to be added
#' @param objectname The object name
#'
#' @export
#' @seealso [pk_add.pk_scales()] for the method for adding `pk_scales` objects
#'   (from [scale_conc()] and [scale_time()]); [pk_add.pk_settings_preprocess()]
#'   for the method for adding `pk_settings_preprocess` objects (from
#'   [settings_preprocess()]); [pk_add.pk_settings_data_info()] for the method
#'   for adding `pk_settings_data_info` objects (from [settings_data_info()]);
#'   [pk_add.pk_settings_optimx()] for the method for adding
#'   `pk_settings_optimx` objects (from [settings_optimx()]);
#'   [pk_add.pk_stat_model()] for the method for adding `pk_stat_model` objects
#'   (from `stat_model()`)
pk_add <- function(pkproto_obj, pk_obj, objectname){
  UseMethod("pk_add")
}

#' Add pkproto object default method
#'
#' @export
pk_add.default <- function(object, pk_obj, objectname){
  stop(paste("No 'pk_add' method exists for object of class",
            object))
}

#'Subtract a `pkproto` object from a `pk` object
#'
#' This is the S3 generic method.
#'
#' @param pkproto_obj The `pkproto` object to be subtracted
#' @param pk_obj The `pk` object to which the `pkproto` object is to be subtracted
#' @param objectname The object name
#'
#' @export
#' @seealso
#'   [pk_subtract.pk_stat_model()] for the method for subtracting `pk_stat_model` objects
#'   (from `stat_model()`)
pk_subtract <- function(pkproto_obj, pk_obj, objectname){
  UseMethod("pk_subtract")
}

#' Subtract pkproto object default method
#'
#' @export
pk_subtract.default <- function(object, pk_obj, objectname){
  stop(paste("No 'pk_subtract' method exists for object of class",
             object))
}

#' Preprocess data generic
#'
#' @export
#' @seealso [preprocess_data.pk()] for the `preprocess_data` method for class [pk()]
preprocess_data <- function(obj, ...){
  UseMethod("preprocess_data", obj)
}

#' Preprocess_data default method
#'
#' @export
preprocess_data.default <- function(obj, ...){
  stop(paste("No 'preprocess_data' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' data_info generic
#'
#' @export
#' @seealso [data_info.pk()] for the `data_info` method for class [pk()]
data_info <- function(obj, ...){
  UseMethod("data_info", obj)
}

#' data_info default method
#'
#' @export
data_info.default <- function(obj, ...){
  stop(paste("No 'data_info' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Prefitting
#'
#' @export
prefit <- function(obj, ...){
  UseMethod("prefit", obj)
}

#' Prefit default method
#'
#' @export
prefit.default <- function(obj, ...){
  stop(paste("No 'prefit' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Fitting
#'
#' This is the S3 generic method for `fit`.
#' @export
#' @seealso [fit.pk()] for the `fit` method for class [pk()]
fit <- function(obj, ...){
  UseMethod("fit", obj)
}

#' Fit default method
#'
#' @export
fit.default <- function(obj, ...){
  stop(paste("No 'fit' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Root mean squared error (RMSE)
#'
#' This is the S3 method generic for `rmse`.
#' @export
#' @seealso [rmse.pk()] for the `rmse` method for class [pk()]
rmse <- function(obj, ...){
  UseMethod("rmse", obj)
}

#' Root mean squared error (RMSE) default method
#'
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
#' @export
#' @seealso [coef_sd.pk()] for the `coef_sd` method for class [pk()]
coef_sd <- function(obj,
                    ...){
  UseMethod("coef_sd", obj)
}

#' Coefficient standard deviation default
#'
#' @param obj An object
#' @export
coef_sd.default <- function(obj, ...){
  stop(paste("No 'coef_sd' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Get status
#'
#' This is the S3 method generic.
#'
#' @param obj An object.
#' @seealso [get_status.pk()] for the `get_status` method for class [pk()]
#' @export
get_status <- function(obj, ...){
  UseMethod("get_status", obj)
}

#' Default method for getting status
#'
#' @export
get_status.default <- function(obj, ...){
  stop(paste("No 'get_status' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Fold error
#'
#' This is the S3 method generic for `fold_errors`.
#' @export
#' @seealso [fold_errors.pk()] for the `fold_errors` method for class [pk()]
fold_errors <- function(obj, ...){
  UseMethod("fold_errors", obj)
}

#' Fold_error default method
#'
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
#' @seealso [check_required_status.pk()] for the method for class [pk()]
#' @export
check_required_status <- function(obj, ...){
  UseMethod("check_required_status", obj)
}

#' Default method for checking required status
#'
#' @export
check_required_status.default <- function(obj, ...){
  stop(paste("No 'check_required_status' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' Get TK stats
#'
#' This is the S3 method generic for get_tkstats(0)
#'
#' @param obj An object.
#' @seealso [get_tkstats.pk()] for the method for class [pk()]
#' @export
get_tkstats <- function(obj, ...){
  UseMethod("get_tkstats", obj)
}

#' Default method for get_tkstats()
#'
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
#' @seealso [compare_models.pk()] for the method for class [pk()]
#' @export
compare_models <- function(obj, ...){
  UseMethod("compare_models", obj)
}

#' Default method for compare_models()
#'
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
#' @seealso [get_data.pk()] for the method for class [pk()]
#' @export
get_data <- function(obj, ...){
  UseMethod("get_data", obj)
}

#' Default method for get_data()
#'
#'@param obj An object
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
#' @seealso [get_nca.pk()] for the method for class [pk()]
#' @export
get_nca <- function(obj, ...){
  UseMethod("get_nca", obj)
}

#' Default method for get_nca()
#'
#'@param obj An object
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
#' @seealso [get_data_info.pk()] for the method for class [pk()]
#' @export
get_data_info <- function(obj, ...){
  UseMethod("get_data_info", obj)
}

#' Default method for get_data_info()
#'
#'@param obj An object
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
#' @seealso [get_prefit.pk()] for the method for class [pk()]
#' @export
get_prefit <- function(obj, ...){
  UseMethod("get_prefit", obj)
}

#' Default method for get_prefit()
#'
#'@param obj An object
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
#' @seealso [get_settings_preprocess.pk()] for the method for class [pk()]
#' @export
get_settings_preprocess <- function(obj, ...){
  UseMethod("get_settings_preprocess", obj)
}

#' Default method for get_settings_preprocess()
#'
#'@param obj An object
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
#' @seealso [get_settings_data_info.pk()] for the method for class [pk()]
#' @export
get_settings_data_info <- function(obj, ...){
  UseMethod("get_settings_data_info", obj)
}

#' Default method for get_settings_data_info()
#'
#'@param obj An object
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
#' @seealso [get_stat_error_model.pk()] for the method for class [pk()]
#' @export
get_stat_error_model <- function(obj, ...){
  UseMethod("get_stat_error_model", obj)
}

#' Default method for get_stat_error_model()
#'
#'@param obj An object
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
#' @seealso [get_data_sigma_group.pk()] for the method for class [pk()]
#' @export
get_data_sigma_group <- function(obj, ...){
  UseMethod("get_data_sigma_group", obj)
}

#' Default method for get_data_sigma_group()
#'
#'@param obj An object
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
#' @seealso [get_settings_optimx.pk()] for the method for class [pk()]
#' @export
get_settings_optimx <- function(obj, ...){
  UseMethod("get_settings_optimx", obj)
}

#' Default method for get_settings_optimx()
#'
#'@param obj An object
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
#' @seealso [get_scale_conc.pk()] for the method for class [pk()]
#' @export
get_scale_conc <- function(obj, ...){
  UseMethod("get_scale_conc", obj)
}

#' Default method for get_scale_conc()
#'
#'@param obj An object
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
#' @seealso [get_scale_time.pk()] for the method for class [pk()]
#' @export
get_scale_time <- function(obj, ...){
  UseMethod("get_scale_time", obj)
}

#' Default method for get_scale_time()
#'
#'@param obj An object
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
#' @seealso [get_error_group.pk()] for the method for class [pk()]
#' @export
get_error_group <- function(obj, ...){
  UseMethod("get_error_group", obj)
}

#' Default method for get_error_group()
#'
#'@param obj An object
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
#' @seealso [get_stat_model.pk()] for the method for class [pk()]
#' @export
get_stat_model <- function(obj, ...){
  UseMethod("get_stat_model", obj)
}

#' Default method for get_stat_model()
#'
#'@param obj An object
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
#' @seealso [get_data_original.pk()] for the method for class [pk()]
#' @export
get_data_original <- function(obj, ...){
  UseMethod("get_data_original", obj)
}

#' Default method for get_data_original()
#'
#'@param obj An object
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
#' @seealso [get_mapping.pk()] for the method for class [pk()]
#' @export
get_mapping <- function(obj, ...){
  UseMethod("get_mapping", obj)
}

#' Default method for get_mapping()
#'
#'@param obj An object
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
#' @seealso [rsq.pk()] for the method for class [pk()]
#' @export
rsq <- function(obj, ...){
  UseMethod("rsq", obj)
}

#' Default method for rsq()
#'
#'@param obj An object
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
#' @seealso [get_winning_model.pk()] for the method for class [pk()]
#' @export
get_winning_model <- function(obj, ...){
  UseMethod("get_winning_model", obj)
}

#' Default method for get_winning_model()
#'
#'@param obj An object
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
#' @seealso [nca.pk()] for the method for class [pk()]
#' @export
nca <- function(obj, ...){
  UseMethod("nca", obj)
}

#' Default method for nca()
#'
#'@param obj An object
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
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
data_summary <- function(obj, ...){
  UseMethod("data_summary", obj)
}

#' Default method for data_summary()
#'
#'@param obj An object
#' @export
data_summary.default <- function(obj, ...){
  stop(paste("No 'data_summary' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_nca()
#'
#' This is the S3 method generic for `get_nca()`
#'
#' `get_nca()` is an alias for `nca`
#'
#' @param obj An object.
#' @seealso [nca.pk()] for the method for class [pk()]
#' @export
get_nca <- function(obj, ...){
  UseMethod("nca", obj)
}

#' Default method for get_nca()
#'
#'@param obj An object
#' @export
get_nca.default <- function(obj, ...){
  stop(paste("No 'get_nca' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}

#' get_data_summary()
#'
#' This is the S3 method generic for get_data_summary()
#'
#' `get_data_summary()` is an alias for `data_summary()`
#'
#' @param obj An object.
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
get_data_summary <- function(obj, ...){
  UseMethod("data_summary", obj)
}

#' Default method for get_data_summary()
#'
#'@param obj An object
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
#' @seealso [data_summary.pk()] for the method for class [pk()]
#' @export
eval_tkstats <- function(obj, ...){
  UseMethod("data_summary", obj)
}

#' Default method for eval_tkstats()
#'
#'@param obj An object
#' @export
eval_tkstats.default <- function(obj, ...){
  stop(paste("No 'eval_tkstats' method exists for object of class",
             paste(class(obj), collapse = ", ")))
}
