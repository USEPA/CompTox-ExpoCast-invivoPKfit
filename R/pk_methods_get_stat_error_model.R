#' Get stat_error_model

#' @param obj A [pk()] object
#' @return A named list of the stat_error_model settings
#' @export
#' @author Caroline Ring
get_stat_error_model.pk <- function(obj){
  out <- obj$stat_error_model
  #convert lists of quosures into "vars(...)"
  out$error_group <- rlang::parse_expr(
    paste0("vars(",
           paste(sapply(out$error_group,
                        function(x) rlang::as_label(x)),
                 collapse = ", "),
           ")")
  )

  out$data_sigma_group <- NULL

  return(out)
}
