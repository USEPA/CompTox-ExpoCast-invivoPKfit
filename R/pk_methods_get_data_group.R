#' Get data grouping
#'
#' @param obj An initialized `pk` object
#' @param ... Additional arguments not currently in use.
#' @return An object of class `call` giving the data grouping as a `dplyr::vars()` specification
#' @export
get_data_group.pk <- function(obj, ...){
  out <- rlang::parse_expr(
    paste0("vars(",
           paste(sapply(obj$data_group,
                        function(x) rlang::as_label(x)),
                 collapse = ", "),
           ")")
  )

  return(out)
}
