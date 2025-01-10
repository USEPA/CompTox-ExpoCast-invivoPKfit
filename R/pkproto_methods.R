#' Is an object pkproto?
#'
#' @param obj An object
#' @return TRUE if `obj` inherits from class `pkproto`; FALSE if not
#' @export
#' @author Caroline Ring
is.pkproto <- function(obj) {
  inherits(obj, "pkproto")
}

#' Is an object class `pk_scales`?
#'
#' @param obj An object
#' @return TRUE if `obj` inherits from class `pk_scales`; FALSE if not
#' @export
#' @author Caroline Ring
is.pk_scales <- function(obj) {
  inherits(obj, "pk_scales")
}
