#' `optimx` optimizer settings
#'
#' @param method The name(s) of optimization methods to be used. See
#'   [optimx::opm()] for options. Default is `"bobyqa"` and `"L-BFGS-B"`.
#' @param hessian Whether to compute the Hessian at the final set of parameters; as in [optimx::opm()].
#' @param control A list of control parameters for the optimizer; see [optimx::opm()] for options and details.
#' @param ... Additional arguments not currently implemented.
#' @return An object of class `pk_settings`.
#' @author Caroline Ring
settings_optimx <- function(method = c("bobyqa", "L-BFGS-B"),
                            hessian = FALSE,
                            control = list(kkt = FALSE, maxit = 1E7),
                            ...) {
  # get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_settings_optimx <- argg
  # set class
  class(this_settings_optimx) <- c(class(this_settings_optimx), "pkproto", "pk_settings_optimx")

  return(this_settings_optimx)
}
