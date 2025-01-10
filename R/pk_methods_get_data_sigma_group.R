#' Get data_sigma_group

#' @param obj A [pk()] object
#' @param newdata Optional: A `data.frame` with new data for which to get the
#'   `data_sigma_group`s. If NULL (the default), then the groups will be
#'   evaluated for the `obj$data`.
#' @param ... Additional arguments. Not currently in use.
#' @return A `factor` vector giving the error SD group ID for each observation,
#'   as the interaction of the factors specified in
#'   `obj$stat_error_model$error_group`.
#' @export
#' @author Caroline Ring
get_data_sigma_group.pk <- function(obj,
                                    newdata = NULL, ...) {
if (is.null(newdata)) {
  newdata <- obj$data
}

  data_sigma_group <- interaction(
    lapply(
      obj$stat_error_model$error_group,
      function(x) {
        rlang::eval_tidy(x, data = newdata)
      }
    )
  )

  return(data_sigma_group)
}
