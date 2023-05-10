#' Add pkproto object
#'
#' @export
pk_add <- function(object, pk_obj, objectname){
  UseMethod("pk_add")
}

#' Add pkproto object default method
#'
#' @export
pk_add.default <- function(object, pk_obj, objectname){
  stop(paste("No 'preprocess_data' method exists for object of class",
            object))
}

#' Preprocess data
#'
#' @export
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
#' @export
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
#' @export
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
#' @param obj An object
#' @export
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

