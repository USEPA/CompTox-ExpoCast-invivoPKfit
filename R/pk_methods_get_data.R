#' Get data
#'
#' Extract pre-processed data from a [pk()] object
#'
#' @param obj A [pk()] object that has been pre-processed
#' @param ... Additional arguments. Currently not in use.
#' @return A `data.frame`: the `data` element of `obj`
#' @export
#' @author Caroline Ring
get_data.pk <- function(obj, ...) {
 # check if data has been preprocessed
  check <- check_required_status(obj = obj,
                                 required_status = status_preprocess)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

return(obj$data)
}
