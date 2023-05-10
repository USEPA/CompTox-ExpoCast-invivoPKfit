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
#'   (from [scale_conc()] and [scale_time()]); [pk_add.pk_data_settings()] for
#'   the method for adding `pk_data_settings` objects (from [settings_data()]);
#'   [pk_add.pk_optimx_settings()] for the method for adding
#'   `pk_optimx_settings` objects (from [settings_optimx()]);
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
preprocess_data <- function(x, ...){
  UseMethod("preprocess_data", x)
}

#' Preprocess_data default method
#'
#' @export
preprocess_data.default <- function(x, ...){
  stop(paste("No 'preprocess_data' method exists for object of class",
             paste(class(x), collapse = ", ")))
}

#' Prefitting
#'
#' @export
prefit <- function(x, ...){
  UseMethod("prefit", x)
}

#' Prefit default method
#'
#' @export
prefit.default <- function(x, ...){
  stop(paste("No 'prefit' method exists for object of class",
             paste(class(x), collapse = ", ")))
}

#' Fitting
#'
#' This is the S3 generic method for `fit`.
#' @export
#' @seealso [fit.pk()] for the `fit` method for class [pk()]
fit <- function(x, ...){
  UseMethod("fit", x)
}

#' Fit default method
#'
#' @export
fit.default <- function(x, ...){
  stop(paste("No 'fit' method exists for object of class",
             paste(class(x), collapse = ", ")))
}

#' Root mean squared error (RMSE)
#'
#' This is the S3 method generic for `rmse`.
#' @export
#' @seealso [rmse.pk()] for the `rmse` method for class [pk()]
rmse <- function(x){
  UseMethod("rmse", x)
}

#' Root mean squared error (RMSE) default method
#'
#' @export
rmse.default <- function(x, ...){
  stop(paste("No 'rmse' method exists for object of class",
             paste(class(x), collapse = ", ")))
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
             paste(class(x), collapse = ", ")))
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
