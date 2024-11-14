#' Get error group

#' @param obj A [pk()] object
#' @param ... Additional arguments. Not in use currently.
#' @return The stat_error_model error grouping
#' @export
#' @author Caroline Ring
get_error_group.pk <- function(obj, ...) {
  out <- rlang::parse_expr(
    paste0("vars(",
           toString(sapply(obj$stat_error_model$error_group,
                           rlang::as_label)
                    ),
           ")")
  )

  return(out)
}
