#' Get data grouping
#'
#' @param obj An initialized `pk` object.
#' @param as_character Logical (Default: `FALSE`). Determines whether to return a
#'  character vector. If set to `FALSE`, a list of expressions containing the `data_group`
#'  variables is returned.
#' @param ... Additional arguments not currently in use.
#' @return An object of class `call` giving the data grouping as a `dplyr::vars()` specification
#' @export
get_data_group.pk <- function(obj, as_character = FALSE, ...) {
  out <- obj$pk_groups$data_group

  if (as_character) {
    out <- vapply(out, rlang::as_string, character(1L))
  }

  return(out)
}
