#' Get settings_optimx

#' @param obj A [pk()] object
#' @param ... Additional arguents not currently in use.
#' @return A named list of the optimx settings
#' @export
#' @author Caroline Ring
get_settings_optimx.pk <- function(obj, ...) {
  out <- obj$settings_optimx

  # convert char vectors into "c(...)
  out$method <- rlang::parse_expr(paste0("c(",
                                         "'",
                                              paste(out$method,
                                                    collapse = "','"),
                                              "')"))

  return(out)
}
