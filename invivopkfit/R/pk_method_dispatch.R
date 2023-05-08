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

#' Preprocess data `pk` method
#'
#' @export
preprocess_data <- function(x){
  UseMethod("preprocess_data", x)
}

#' Preprocess_data default method
#'
#' @export
preprocess_data.default <- function(x){
  stop(paste("No 'preprocess_data' method exists for object of class",
             paste(class(x), collapse = ", ")))
}
