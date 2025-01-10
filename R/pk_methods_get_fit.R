#' Get fits from a `pk` object
#'
#' Get the [optimx::optimx()] output from a fitted `pk` object
#'
#' This function returns the object(s) returned by [optimx::optimx()] for the
#' specified model(s) and method(s), for a fitted `pk` object. See
#' [optimx::optimx()] for details. Briefly, an `optimx` object is a `data.frame`
#' with one row for each method used, and variables that give the optimized
#' values for each parameter, along with several diagnostic variables (e.g. the
#' objective function value at the optimized parameter values; the number of
#' function evaluations/iterations; an integer code describing convergence
#' status). The object will have attributes `details` (providing
#' any messages returned by the methods) and `npar` (the number of parameters
#' optimized).
#'
#' @param obj A [pk] object.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions. If NULL (the default), predictions will be returned for
#'   all of the models in `obj$stat_model`.
#' @param ... Additional arguments. Not in use.
#' @return A named list of objects of class `optimx`, named for the models in
#'   `model`. As described in [optimx::optimx()]  If only one model is
#'   specified, the return value will still be a list, but with only one
#'   element.
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @export
#' @family methods for fitted pk objects
get_fit.pk <- function(obj,
                       model = NULL,
                       ...) {

  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model)) model <- names(obj$stat_model)

  model_ok <- check_model(obj = obj, model = model)

  tmp <- subset(obj$fit, model %in% model)

  return(tmp)
}
