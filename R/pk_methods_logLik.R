#' Log-likelihood
#'
#' Extract log-likelihood(s) from a fitted `pk` object
#'
#' For details on how the log-likelihood is calculated, see [log_likelihood()].
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then log-likelihoods will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `Conc_SD`, `N_Subjects`. Before log-likelihood is calculated,
#'   `Time` will be transformed according to the transformation in
#'   `obj$scales$time` and `Conc` will be transformed according to the
#'   transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate log-likelihood. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to calculate log-likelihoods. If NULL (the default),
#'   log-likelihoods will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @return A named list of numeric vectors. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a numeric vector with as many
#'   elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the log-likelihood of the model fitted by the
#'   corresponding method, calculated for the data in `newdata` if any has been
#'   specified, or for the data in `obj$data` if `newdata` is `NULL`.
#' @export
#' @importFrom stats logLik
#' @author Caroline Ring
logLik.pk <- function(obj,
                      newdata = NULL,
                      model = NULL,
                      method = NULL,
                      negative = FALSE,
                      force_finite = FALSE){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

  #transform time needed
  #first, default to identity transformation if none is specified
  if(is.null(obj$scales$time$new_units)){
    obj$scales$time$new_units <- "identity"
  }

  from_units <- unique(newdata$Time.Units)
  to_units <- ifelse(obj$scales$time$new_units %in% "identity",
                     from_units,
                     obj$scales$time$new_units)


  if(obj$scales$time$new_units %in% "auto"){
    to_units <- auto_units(y = newdata$Time,
                           from = from_units)
  }

  newdata$Time_trans <- tryCatch(convert_time(x = newdata$Time,
                                              from = from_units,
                                              to = to_units,
                                              inverse = FALSE),
                                 error = function(err){
                                   warning(paste("invivopkfit::logLik.pk():",
                                                 "Error in transforming time in `newdata` using convert_time():",
                                                 err$message))
                                   return(NA_real_)
                                 })

  #Apply concentration transformation
  newdata$Conc_trans <- rlang::eval_tidy(obj$scales$conc$expr,
                                         data = cbind(newdata,
                                                      data.frame(".conc" = newdata$Conc)))
  newdata$Conc_SD_trans <- rlang::eval_tidy(obj$scales$conc$expr,
                                            data = cbind(newdata,
                                                         data.frame(".conc" = newdata$Conc_SD)))

  sapply(model, function(this_model){
    #if there is newdata, then we have to evaluate the log-likelihood
    #get model parameters
    coefs <- coef(obj = obj,
                  model = this_model,
                  method = method)[[1]] #we have to put [[1]] because for one model,coef.pk() returns a one-element list
    #for each row of model parameters, evaluate log-likelihood
    ll <- apply(coefs,
                1,
                function(this_coef_row){
                  log_likelihood(par = as.list(this_coef_row),
                                 fitdata = newdata,
                                 data_sigma_group = obj$stat_error_model$data_sigma_group,
                                 modelfun = obj$stat_model[[this_model]]$conc_fun,
                                 scales_conc = obj$scales$conc,
                                 negative = negative,
                                 force_finite = force_finite)
                })
    #set attribute "df", the number of parameters optimized for this model
    attr(ll, which = "df") <- attr(obj$stat_model[[this_model]]$fit, "npar")
    #set attributes "nobs", the number of observations in `newdata`
    attr(ll, which = "nobs") <- nrow(newdata)
    ll
  },
  simplify = FALSE,
  USE.NAMES = TRUE)

}
