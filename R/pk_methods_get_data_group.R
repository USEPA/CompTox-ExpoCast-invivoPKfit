#' Get data grouping
#'
#' @param obj An initialized `pk` object
#' @param ... Additional arguments not currently in use.
#' @return An object of class `call` giving the data grouping as a `dplyr::vars()` specification
#' @export
get_data_group.pk <- function(obj, ...) {
  out <- rlang::parse_expr(
    paste0("vars(",
           paste(
             toString(
               sapply(obj$data_group,
                      rlang::as_label)
             ),
           ")"
           )
    )
  )
  return(out)
}
