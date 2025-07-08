#' Get settings_data_info

#' @param obj A [pk()] object
#' @param ... Additional arguments. Currently not implemented.
#' @return A named list of the data_info settings
#' @export
#' @author Caroline Ring
get_settings_data_info.pk <- function(obj, ...) {
  out <- obj$settings_data_info
  # convert lists of quosures into "vars(...)"
  out$summary_group <- rlang::parse_expr(
    paste0("vars(",
           toString(sapply(out$summary_group, rlang::as_label)),
           ")")
  )

  return(out)
}
