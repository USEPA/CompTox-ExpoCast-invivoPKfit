#' Akaike information criterion
#'
#' Get the Akaike information criterion (AIC) for a fitted `pk` object
#'
#' The AIC is calculated from the log-likelihood (LL) as follows:
#' \deqn{\textrm{AIC} = -2\textrm{LL} + k n_{par}}
#'
#' where \eqn{n_{par}} is the number of parameters in the fitted model, and
#' \eqn{k = 2} for the standard AIC.
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then log-likelihoods will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Time.Units`, `Dose`,
#'   `Route`,`Media`, `Conc`, `Detect`, `N_Subjects`. Before log-likelihood is
#'   calculated, `Time` will be transformed according to the transformation in
#'   `obj$scales$time` and `Conc` will be transformed according to the
#'   transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate log-likelihood. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate AICs. If NULL (the default),
#'   log-likelihoods will be returned for all of the models in
#'   `obj$settings_optimx$method`.
#' @param exclude Logical: `TRUE` to compute the AIC after removing any
#'   observations in the data marked for exclusion (if there is a variable
#'   `exclude` in the data, an observation is marked for exclusion when `exclude
#'   %in% TRUE`). `FALSE` to include all observations, regardless of exclusion
#'   status. Default `TRUE`.
#' @param k Default 2. The `k` parameter in the log-likelihood formula (see
#'   Details).
#' @return A named list of numeric vectors. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a numeric vector with as many
#'   elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the AIC of the model fitted by the corresponding method,
#'   using the data in `newdata`.
#' @family fit evaluation metrics
#' @family log likelihood functions
#' @family methods for fitted pk objects
#' @importFrom stats AIC
#' @export
#' @author Caroline Ring
AIC.pk <- function(obj,
                   newdata = NULL,
                   model = NULL,
                   method = NULL,
                   exclude = TRUE,
                   k = 2){
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = 5)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method
  if(is.null(newdata)) newdata <- obj$data
  #get log-likelihoods
  ll <- logLik(obj = obj,
               newdata = newdata,
               model = model,
               method = method,
               negative = FALSE,
               force_finite = FALSE,
               exclude = exclude)
  #get number of parameters (excluding any constant, non-optimized parameters)
  AIC <- sapply(model,
                function(this_model){
                  npar <- attr(ll[[this_model]], "df")
                  -2*ll[[this_model]] + k* npar
                },
                simplify = FALSE,
                USE.NAMES = TRUE)
  return(AIC)
}
