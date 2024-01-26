#' Get scale_time
#'
#' Extract time scale/transformation instructions from a [pk()] object
#'
#' @param obj A [pk()] object
#' @param ... Additional arguments not in use.
#' @return A `list`: `obj$scales$time`
#' @author Caroline Ring
#' @export
get_scale_time.pk <- function(obj, ...){
  return(obj$scales$time)
}
