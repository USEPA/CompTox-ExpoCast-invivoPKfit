#' Get prefit
#'
#' Extract pre-fitting results from a [pk()] object
#'
#' @param obj A [pk()] object that has had `do_prefit()` run on it
#' @param ... Additional arguments. Currently not in use.
#' @return A list of `data.frame`s in the object's `prefit` element.
#' @export
#' @author Caroline Ring
get_prefit.pk <- function(obj,
                          ...) {

  # check if data has been prefit
  check <- check_required_status(obj = obj,
                                 required_status = status_prefit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  return(obj$prefit)
}
