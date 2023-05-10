#' Root mean squared error
#'
#' Extract root mean squared error of a fitted `pk` object
#'
#' RMSE is calculated using the following formula, to properly handle summary
#' data:
#'
#' \deqn{
#' \sqrt{
#' \frac{1}{N}
#'  \sum_{i=1}^G \left( (n_i - 1) s_i^2 +
#'   n_i ^2 - 2 n_i \bar{y}_i \mu_i + \mu_i^2 \right)
#'    }
#' }
#'
#' In this formula, there are \eqn{G} observations, each of which may be for one
#' subject or for multiple subjects.
#'
#' - \eqn{n_i} is the number of subjects for observation \eqn{i}.
#' - \eqn{\bar{y}_i} is the sample mean concentration for observation \eqn{i}.
#' - \eqn{s_i} is the sample standard deviation of concentrations for observation \eqn{i}.
#' - \eqn{\mu_i} is the model-predicted concentration for observation \eqn{i}.
#'
#' \eqn{N} is the grand total of subjects across observations:
#'
#' \deqn{N = \sum_{i=1}^G n_i}
#'
#' For the non-summary case (\eqn{N} single-subject observations, with all
#' \eqn{n_i = 1}, \eqn{s_i = 0}, and \eqn{\bar{y}_i = y_i}), this formula
#' reduces to the familiar RMSE formula
#'
#' \deqn{\sqrt{\frac{1}{N} \sum_{i=1}^N (y_i - \mu_i)^2}}
#'
#' # Left-censored data
#'
#' If the observed value is censored, and the predicted value is less than the
#' reported LOQ, then the predicted value is (temporarily) set equal to the LOQ,
#' for an effective error of zero.
#'
#' If the observed value is censored, and the predicted value is greater than
#' the reported LOQ, the the observed value is treated as the reported LOQ (so
#' that the effective error is the difference between the LOQ and the predicted
#' value).
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute RMSsE. If NULL (the default), then RMSEs will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`, and `Media`. `Time`
#'   will be transformed according to the transformation in `obj$scales$time`
#'   before RMSEs are calculated.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate RMSEs. If NULL (the default), RMSEs will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate RMSEs. If NULL (the default),
#'   RMSEs will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC). Currently, only `type = "conc"` is
#'   implemented.
#' @return A named list of numeric vectors. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a numeric vector with as many
#'   elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the root mean squared error of the model fitted by the
#'   corresponding method, using the data in `newdata`. These RMSEs are
#'   concentrations in the same units as `obj$data$Conc.Units`; any
#'   concentration transformations (in `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
rmse.pk <- function(obj,
                    newdata = NULL,
                    model = NULL,
                    method = NULL){

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

  #get predicted concentrations
  preds <- predict(obj,
                      newdata = newdata,
                      model = model,
                      method = method,
                      type = "conc")

  n_subj <- newdata$N_Subjects
  obs <- newdata$Conc
 obs_sd <- newdata$Conc_SD
  n_tot <- sum(n_subj)


  sapply(preds,
         function(this_pred){
           apply(this_pred,
                 2,
                 function(x) {
                   x <- ifelse(newdata$Detect %in% FALSE &
                                 x < obs,
                               obs, #will be LOQ for detect == FALSE
                               x)

                   sqrt(
                     sum(
                       (n_subj -1) * obs_sd^2 +
                         n_subj * obs^2 -
                         2 * n_subj*obs*x +
                         x^2
                     )/n_tot
                   )
                 }
           )
         },
         simplify = FALSE,
         USE.NAMES = TRUE)
}
