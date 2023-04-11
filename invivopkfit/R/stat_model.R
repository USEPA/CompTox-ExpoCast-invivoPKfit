#' Add PK model(s) to be fitted
#'
#' @param obj A `pk` object
#' @param model The name(s) of models to be fitted. Options are `"flat"`,
#'   `"1comp"`, and `"2comp"`. Other models are not currently implemented.
#' @param method The name(s) of optimization methods to be used. See
#'   [optimx::optimx()] for options. Default is `"bobyqa"`.
#' @param itnmax The maximum number of iterations; as in [optimx::optimx()].
#' @param hessian Whether to compute the Hessian at the final set of parameters; as in [optimx::optimx()].
#' @param control A list of control parameters for the optimizer; see [optimx::optimx()] for options and details.
stat_model.pk <- function(obj,
                       model = c("flat", "1comp", "2comp"),
                       method = "bobyqa",
                       itnmax = 1e6,
                       hessian = FALSE,
                       control = list(kkt = FALSE)){
  #add model fitting options for each model to be fitted
  for (this_model in model){
    obj$models[[this_model]][["settings"]] <- list(method = method,
                                     itnmax = itnmax,
                                     hessian = hessian,
                                     control = control)
  }

  return(obj)
}
