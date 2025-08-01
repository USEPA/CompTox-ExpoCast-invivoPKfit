#' Get error group

#' @param obj A [pk()] object.
#' @param as_character Logical (Default: `FALSE`). Determines whether to return a
#'  character vector. If set to `FALSE`, a list of expressions containing the `data_group`
#'  variables is returned.
#' @param ... Additional arguments. Not in use currently.
#' @return The stat_error_model error grouping
#' @export
#' @author Caroline Ring
get_error_group.pk <- function(obj, as_character = FALSE, ...) {
  out <- obj$pk_groups$error_group

  if (as_character) {
    out <- vapply(out, rlang::as_string, character(1L))
  }


  return(out)
}
