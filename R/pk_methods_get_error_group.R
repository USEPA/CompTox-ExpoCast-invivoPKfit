#' Get error group

#' @param obj A [pk()] object
#' @return The stat_error_model error grouping
#' @export
#' @author Caroline Ring
get_error_group.pk <- function(obj){
  out <- rlang::parse_expr(
    paste0("vars(",
           paste(sapply(obj$stat_error_model$error_group,
                        function(x) rlang::as_label(x)),
                 collapse = ", "),
           ")")
  )

  return(out)
}
