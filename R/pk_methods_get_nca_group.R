#' Get nca_group

#' @param obj A [pk()] object.
#' @param as_character Logical (Default: `FALSE`). Determines whether to return a
#'  character vector. If set to `FALSE`, a list of expressions containing the `data_group`
#'  variables is returned.
#' @param ... Additional arguments. Currently not implemented.
#' @return A named list of the data_info settings
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
get_nca_group.pk <- function(obj, as_character = FALSE, ...) {
  out <- obj$pk_groups$nca_group

  if (as_character) {
    out <- vapply(out, rlang::as_string, character(1L))
  }

  return(out)
}
