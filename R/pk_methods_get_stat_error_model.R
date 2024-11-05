#' Get stat_error_model

#' @param obj A [pk()] object
#' @param ... Additional arguments.
#' @return A named list of the stat_error_model settings
#' @export
#' @author Caroline Ring
get_stat_error_model.pk <- function(obj, ...){
  out <- obj$stat_error_model
  #convert lists of quosures into "vars(...)"
  out$error_group <- rlang::parse_expr(
    paste0("vars(",
           paste(toString(out$error_group,
                          rlang::as_label),
           ")"
           )
    )
  )

  out$data_sigma_group <- NULL

  return(out)
}
