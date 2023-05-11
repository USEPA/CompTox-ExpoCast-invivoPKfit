#' Bayesian information criterion
#'
#' Get the Bayesian information criterion (AIC) for a fitted `pk` object
#'
#' The BIC is calculated from the log-likelihood (LL) as follows:
#' \deqn{\textrm{BIC} = -2\textrm{LL} + \log(n_{obs}) n_{par}}
#'
#' where \eqn{n_{par}} is the number of parameters in the fitted model.
#'
#' Note that the BIC is just the AIC with \eqn{k = \(n_{obs})}.
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then BICs will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `N_Subjects`. Before log-likelihood is calculated, `Time` will be
#'   transformed according to the transformation in `obj$scales$time` and `Conc`
#'   will be transformed according to the transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate BIC. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to calculate BICs. If NULL (the default),
#'   log-likelihoods will be returned for all of the methods in
#'   `obj$optimx_settings$method`.
#' @return A named list of numeric vectors. There is one list element
#'   named for each model in `obj`'s [stat_model()] element, i.e. each PK model
#'   that was fitted to the data. Each list element is a numeric vector with as
#'   many elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the BIC of the model fitted by the corresponding method,
#'   using the data in `newdata`.
#' @seealso [AIC.pk(0)], [logLik.pk()]
#' @importFrom stats BIC
#' @export
#' @author Caroline Ring
BIC.pk <- function(obj,
                   newdata = NULL,
                   model = NULL,
                   method = NULL){
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = 4)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

  BIC <- AIC(obj = obj,
             newdata = newdata,
             model = model,
             method = method,
             k = log(nrow(newdata)))
  return(BIC)
}
