#' Get mapping
#'
#' @param obj A [pk()] object
#' @return A list of `quosure`s -- the `mapping` element of `obj`
#' @export
#' @author Caroline Ring
get_mapping.pk <- function(obj){
  obj$mapping
}
