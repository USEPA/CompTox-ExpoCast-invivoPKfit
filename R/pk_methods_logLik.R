#' Log-likelihood
#'
#' Extract log-likelihood(s) from a fitted `pk` object
#'
#' For details on how the log-likelihood is calculated, see [log_likelihood()].
#'
#' # New levels in `newdata`
#'
#' The log-likelihood requires an error variance for each observation. Depending
#' on the error model used for the model fit, each observation may have a
#' different corresponding error variance, based on its unique combinations of
#' levels of factor variables as specified in [stat_error_model()] when setting
#' up the [pk()] object. (The error model specified in the `pk` object can be
#' viewed using [get_error_group()]).
#'
#' If you are supplying new data in `newdata`, then the unique combinations of
#' the factor levels for the new observations will be used to find the matching
#' error hyperparameter. If the new data contains new combinations of factor
#' levels not found in the data used to fit the model, then the following
#' procedure will be used to calculate the log-likelihood for the observations
#' with new levels:
#'
#' If there are \eqn{i = 1, 2, ... G} groups with unique combinations of the
#' factor levels in the original data, then there are corresponding error
#' standard deviations \eqn{\sigma_i = \sigma_1, \sigma_2, ..., \sigma_G}.
#'
#' Each observation with a new combination of factor levels will have \eqn{G}
#' differerent log-likelihoods computed, as though it were part of each of the
#' \eqn{G} existing groups. Then, the average of these \eqn{G} log-likelihoods
#' will be taken and assigned to the observation. In effect, each observation
#' with a new level is treated as though it is equally likely to belong to any
#' of the existing groups.
#'
#' # Scaling and transformation of concentration variables in `newdata`
#'
#' This function differs from some of the other methods for a fitted [pk()]
#' object that accept `newdata`, in that there is no `use_scale_conc` argument
#' for [logLik.pk()]. You cannot specify a different scaling/transformation for
#' concentration variables -- you have to use the same scaling/transformation
#' that was used to fit the model. This is because the log-likelihood depends on
#' the fitted values of the error variance hyperparameters, and those are valid
#' only for the transformation of the concentration data that was used to fit
#' the model.
#'
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then log-likelihoods will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Time.Units`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `Conc_SD`, `N_Subjects`. `Time` will be transformed according to
#'   the transformation in `obj$scales$time` before making predictions. `Conc`
#'   will be transformed according to the transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate log-likelihood. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to calculate log-likelihoods. If NULL (the default),
#'   log-likelihoods will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param exclude Logical: `TRUE` to compute the log-likelihood excluding any
#'   observations in the data marked for exclusion (if there is a variable
#'   `exclude` in the data, an observation is marked for exclusion when `exclude
#'   %in% TRUE`). `FALSE` to include all observations when calculating the
#'   log-likelihood, regardless of exclusion status. Default `TRUE`.
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
#' @family fit evaluation metrics
#' @family log likelihood functions
#' @family methods for fitted pk objects
logLik.pk <- function(obj,
                      newdata = NULL,
                      model = NULL,
                      method = NULL,
                      negative = FALSE,
                      force_finite = FALSE,
                      exclude = TRUE){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data



  #check variables in newdata
#including error grouping variables
  err_grp_vars <- sapply(eval(get_error_group(obj)),
                         rlang::as_label)

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = union(c("Time",
                                           "Dose",
                                           "Route",
                                           "Media",
                                           "Conc",
                                           "Detect",
                                           "Conc_SD",
                                           "N_Subjects"),
                                err_grp_vars),
                              exclude = exclude)

  if(exclude %in% TRUE){
    if("exclude" %in% names(newdata)){
    newdata <- subset(newdata, exclude %in% FALSE)
    }
  }else{
    newdata <- newdata
  }

  #scale time if needed
  if(!("Time_trans" %in% names(newdata))){
    newdata$Time_trans <- convert_time(x = newdata$Time,
                                       from = newdata$Time.Units,
                                       to = obj$scales$time$new_units)
  }

  #get transformations to apply
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = TRUE)

  #remove any excluded observations & corresponding predictions, if so specified
  if(exclude %in% TRUE){
    if("exclude" %in% names(newdata)){
      newdata <- subset(newdata,
                        exclude %in% FALSE)

    }
  }

  #evalaute log likelihoods for each model using sapply()
  sapply(model, function(this_model){
    #get model parameters
    coefs <- coef(obj = obj,
                  model = this_model,
                  method = method)[[1]] #we have to put [[1]] because for one model,coef.pk() returns a one-element list
    #for each row of model parameters, evaluate log-likelihood
    ll <- apply(coefs,
                1,
                function(this_coef_row){

                  data_sigma_group <- as.character(interaction(
                    lapply(
                      obj$stat_error_model$error_group,
                      function(x){
                        rlang::eval_tidy(x, data = newdata)
                      }
                    )
                  ))

                  #apply levels from original data
                  data_sigma_group <- factor(data_sigma_group,
                                             levels = levels(obj$prefit$stat_error_model$data_sigma_group))

                 #data from new levels will be treated as equally likely to come from any of the existing levels
                  log_likelihood(par = as.list(this_coef_row),
                                 data = newdata,
                                 data_sigma_group = data_sigma_group,
                                 modelfun = obj$stat_model[[this_model]]$conc_fun,
                                 dose_norm = conc_scale$dose_norm,
                                 log10_trans = conc_scale$log10_trans,
                                 negative = negative,
                                 force_finite = force_finite)
                })
    #set attribute "df", the number of parameters optimized for this model
    attr(ll, which = "df") <- attr(obj$fit[[this_model]], "npar")
    #set attributes "nobs", the number of observations in `newdata`
    attr(ll, which = "nobs") <- nrow(newdata)
    ll
  },
  simplify = FALSE,
  USE.NAMES = TRUE)

}
