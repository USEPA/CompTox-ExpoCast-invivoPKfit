#' Check whether an object is of class `pk`
#'
#' @param obj The object whose class is to be tested
#' @return TRUE if the object inherits from class `pk`, FALSE if it does not
#' @export
#' @author Caroline Ring
is.pk <- function(obj){
  return(inherits(obj, "pk"))
}

#' Check whether an object is of class `pk_faceted`
#'
#' @param obj The object whose class is to be tested
#' @return TRUE if the object inherits from class `pk`, FALSE if it does not
#' @export
#' @author Caroline Ring
is.pk_faceted <- function(obj){
  return(inherits(obj, "pk_faceted"))
}
