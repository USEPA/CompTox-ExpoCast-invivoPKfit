#' Get settings_preprocess

#' @param obj A [pk()] object
#' @return A named list of the preprocessing settings
#' @export
#' @author Caroline Ring
get_settings_preprocess.pk <- function(obj){
  out <- obj$settings_preprocess
  #convert char vectors into "c(...)
  out$routes_keep <- rlang::parse_expr(paste0("c(",
                                              "'",
                                              paste(out$routes_keep,
                                                    collapse = "','"),
                                              "')"))

  out$media_keep <- rlang::parse_expr(paste0("c(",
                                             "'",
                                              paste(out$media_keep,
                                                    collapse = "','"),
                                              "')"))

  #convert lists of quosures into "vars(...)"
  out$loq_group <- rlang::parse_expr(
    paste0("vars(",
                         paste(sapply(out$loq_group,
                          function(x) rlang::as_label(x)),
                          collapse = ", "),
                         ")")
  )

  out$sd_group <- rlang::parse_expr(
    paste0("vars(",
                          paste(sapply(out$sd_group,
                                       function(x) rlang::as_label(x)),
                                collapse = ", "),
                          ")")
  )



  return(out)
}
