#' Get data groupingh
#'
#' @param obj An initialized `pk` object
#' @param ... Additional arguments not currently in use.
#' @return An object of class `call` giving the data grouping as a `dplyr::vars()` specification
#' @export
get_data_group.pk <- function(obj, ...){
  obj$data_group
}
