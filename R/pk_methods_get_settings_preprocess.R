#' Get settings_preprocess

#' @param obj A [pk()] object
#' @param ... Additional arguments. Currently not in use.
#' @return A named list of the preprocessing settings
#' @export
#' @author Caroline Ring
get_settings_preprocess.pk <- function(obj, ...){
  out <- obj$settings_preprocess
  #convert char vectors into "c(...)
  out$routes_keep <- rlang::parse_expr(paste0("c(",
                                              "'",
                                              toString(out$routes_keep),
                                              "')"))

  out$media_keep <- rlang::parse_expr(paste0("c(",
                                             "'",
                                              toString(out$media_keep),
                                              "')"))

  #convert lists of quosures into "vars(...)"
  out$loq_group <- rlang::parse_expr(
    paste0("vars(",
           toString(sapply(out$loq_group,
                           function(x) rlang::as_label(x))),
           ")")
  )

  out$sd_group <- rlang::parse_expr(
    paste0("vars(",
           toString(sapply(out$sd_group, function(x) rlang::as_label(x))),
           ")")
  )



  return(out)
}
