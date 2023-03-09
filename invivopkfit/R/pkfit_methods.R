#' Create a `pkfit` object
#'
#'
pkfit <- function(data = NULL,
                  data_info = )


#' Create a `pkfit` object from data
#'
#' Given a single-row data.frame of fit outputs and a data.frame of data,
#' produce an object of class `pkfit`
#'
#' @param fit A single-row `data.frame`. One row of the `data.frame` produced
#'   e.g. by [fit_all()].
#' @param dat A `data.frame` of concentration vs. time data as passed to
#'   [fit_all()], as produced by [preprocess_data()].
#' @return An pbject of class `pkfit`. This is a `list` with named elements
#'   `data`, `data_info`, `data_trans`, `model_info`, `fit_params`, `fit_gof`, `fit_info`,
#'   `fit_control`.
#' @export
#' @author Caroline Ring
as.pkfit <- function(fit, dat){
  fit <- as.data.frame(fit)
  dat <- as.data.frame(dat)

  #get the analysis type
  pkfit_analysis_type <- unique(fit$Analysis_Type)

  #get the data info
  pkfit_dat_info <- as.list(unique(fit[c("DTXSID",
                                  "Species",
                                  "Studies.Analyzed")]))
  #get a list of studies analyzed
  pkfit_dat_info$Studies_Analyzed <- strsplit(unique(fit$Studies.Analyzed),
                                              split = ", ")[[1]]
  #get a list of routes analyzed
  tmp <- strsplit(
    strsplit(unique(fit$N_Routes),
                                  split = " ; ")[[1]],
    split = ": ")
  tmp_routes <- sapply(tmp, function(x) x[1])
  tmp_n <- sapply(tmp, function(x) as.numeric(x[2]))
  pkfit_dat_info$Routes_Analyzed <- tmp_routes[tmp_n > 0]

  #pull the appropriate subset of data
  pkfit_dat <- subset(dat,
                      DTXSID %in% pkfit_dat_info$DTXSID &
                      Species %in% pkfit_dat_info$Species &
                      Route %in% pkfit_dat_info$Routes_Analyzed &
                        Study %in% pkfit_dat_info$Studies_Analyzed)
  if(!("Conc_SD" %in% names(pkfit_dat)) &
     "Value_SD" %in% names(pkfit_dat)){
    pkfit_dat$Conc_SD <- pkfit_dat$Value_SD
    pkfit_dat$Conc_SD_Dose <- pkfit_dat$Conc_SD/pkfit_dat$Dose
  }

  #grab the data transformation info
  pkfit_trans <- unique(fit[c("fit_conc_dose",
                               "fit_log_conc",
                               "rescale_time",
                               "time_units_fitted")])

  #grab the fitted model info
  pkfit_model <- as.list(unique(fit[c("model",
                               "model.type")]))

  #grab the fitted parameter info
  pkfit_params <- fit[c("param_name",
                         "param_units",
                         "optimize_param",
                         "use_param",
                         "lower_bound",
                         "lower_bound_msg",
                         "upper_bound",
                         "upper_bound_msg",
                         "start_value",
                         "start_value_msg",
                         "Fitted mean",
                         "Fitted std dev")]

  #grab the goodness of fit metrics (AIC and log-likelihood)
  pkfit_gof <- unique(fit[c("AIC",
                             "LogLikelihood")])

  #grab the optimizer arguments
  pkfit_optim_args <- list("method" = unique(fit$method))
  if(!("itnmax" %in% names(fit))){
    pkfit_optim_args$itnmax <- 1e6
  }

  pkfit_optim_out <- unique(fit[c("message",
                               "fevals",
                               "convcode",
                               "niter")])
  #grab the optimizer control info
  pkfit_optim_control <- unique(fit[grep(x = names(fit),
                                          pattern = "control_",
                                          fixed = TRUE,
                                          value = TRUE)])
  names(pkfit_optim_control) <- gsub(x = names(pkfit_optim_control),
                                     pattern = "control_",
                                     replacement = "",
                                     fixed = TRUE)
  pkfit_optim_control <- as.list(pkfit_optim_control)



  # construct the pkfit object
  obj <- list("data" = pkfit_dat,
              "data_info" = pkfit_dat_info,
              "data_trans" = pkfit_trans,
              "analysis_type" = pkfit_analysis_type,
              "model_info" = pkfit_model,
              "fit_params" = pkfit_params,
              "fit_gof" = pkfit_gof,
              "fit_args" = pkfit_optim_args,
              "fit_control" = pkfit_optim_control,
              "fit_info" = pkfit_optim_out
             )

  class(obj) <- c(class(obj), "pkfit")

  return(obj)
}

#' Check whether an object is of class pkfit
#'
#' @param obj The object whose class is to be tested
#' @return TRUE if the object conforms to expectations for class
#'   `pkfit`, FALSE if it does not
#' @export
#' @author Caroline Ring
is.pkfit <- function(obj){
  expected_names <- c("data",
                     "data_info",
                     "data_trans",
                     "analysis_type",
                     "model_info",
                     "fit_params",
                     "fit_gof",
                     "fit_args",
                     "fit_control",
                     "fit_info")
  expected_classes <- c("data" = "data.frame",
                        "data_info" = "list",
                        "data_trans" = "data.frame",
                        "analysis_type" = "character",
                        "model_info" = "list",
                        "fit_params" = "data.frame",
                        "fit_gof" = "data.frame",
                        "fit_args" = "list",
                        "fit_control"= "list",
                        "fit_info" = "data.frame")

  check <- inherits(obj, "pkfit") &
    all(expected_names %in% names(obj))

  #if we're OK so far, check whether each named element inherits from the
  #expected class
  if(check %in% TRUE){
   check <- check &
     all(
       sapply(expected_names,
          function(x){
            inherits(obj[[x]],
                     expected_classes[x])
          })
     )
  }

  return(check)
}

#' Extract fitted parameter values as a list
#'
#' This function differs from [coefficients.pkfit()] only in that it returns a
#' named *list* of fitted model coefficients, rather than a named *vector*. This
#' difference arises because the model functions (e.g. [cp_flat()],
#' [cp_1comp()], [cp_2comp()]) require a named *list* of parameters, so it is
#' convenient to have a method to access the fitted coefficients in list form.
#'
#' @param obj A `pkfit` object
#' @return A named list of the fitted parameter values
#' @export
#' @author Caroline Ring
get_params.pkfit <- function(obj){
  #get the params
  params <- obj$fit_params[["Fitted mean"]]
  names(params) <- obj$fit_params[["param_name"]]
  params <- as.list(params)

  #keep only finite param values
  params <- params[sapply(params,
                          is.finite)]

  #convert time constants back to hours if necessary
  if(!(obj$data_trans$time_units_fitted %in% "hours")){
    time_const <- intersect(names(params),
                            c("kelim",
                              "kgutabs",
                              "k12",
                              "k21"))

  for(this_tau in time_const){
    params[[this_tau]] <- convert_time(params[[this_tau]],
                                   from = obj$data_trans$time_units_fitted,
                                   to = "hours",
                                   inverse = TRUE)
  }
  }
  return(params)
}

#' Get coefficients of a `pkfit` object
#'
#' @param obj A fitted `pkfit` object
#' @return A named numeric vector of the fitted model coefficients
#' @export
#' @author Caroline Ring
coef.pkfit <- function(obj){
  params <- get_params.pkfit(obj)
  return(unlist(params))
}

#' Get coefficients of a `pkfit` object
#'
#' Alias for [coef.pkfit()]
#'
#' @param obj A fitted `pkfit` object
#' @return A named numeric vector of the fitted model coefficients
#' @export
#' @author Caroline Ring
coefficients.pkfit <- function(obj){
  coef.pkfit(obj)
}

#' Extract fitted parameter standard deviations as a list
#'
#' @param obj A `pkfit` object
#' @return A named list of the fitted parameter standard deviations. The names
#'   are the parameter names with "_sd" appended.
#' @export
#' @author Caroline Ring
get_params_sd.pkfit <- function(obj){
  #get the params
  params <- obj$fit_params[["Fitted std dev"]]
  names(params) <- obj$fit_params[["param_name"]]
  params <- as.list(params)

  #keep only finite param values
  params <- params[sapply(params,
                          is.finite)]

  #convert time constants back to hours if necessary
  if(!(obj$data_trans$time_units_fitted %in% "hours")){
    time_const <- intersect(names(params),
                            c("kelim",
                              "kgutabs",
                              "k12",
                              "k21"))

  for(this_tau in time_const){
    params[[this_tau]] <- convert_time(params[[this_tau]],
                                       from = obj$data_trans$time_units_fitted,
                                       to = "hours",
                                       inverse = TRUE)
  }
  }
if(length(params)>0){
  names(params) <- paste(names(params),
                         "sd",
                         sep = "_")
}
  return(params)
}

#' Get predicted concentrations for a `pkfit` object

#' @param obj The `pkfit` object for which to return predicted concentrations
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), and `Media` ("blood" or "plasma").
#' @param what Character: "conc" to predict concentrations, "auc" to predict AUC
#'   (area under the concentration-time curve). Not case sensitive.
#' @return A numeric vector of predicted concentrations with length equal to
#'   `nrows(newdata)`
#' @export
#' @author Caroline Ring
predict.pkfit <- function(obj,
                          newdata = NULL,
                          what = "conc"){
  #get the model function to be evaluated
  if(tolower(what) %in% "conc"){
    modelfun <- get_model_function(model = obj$model_info$model,
                                   model.type = obj$model$model.type)
  }else if(tolower(what) %in% "auc"){
    modelfun <- get_model_auc_function(model = obj$model_info$model,
                                   model.type = obj$model$model.type)
  }
    #get the params
   params <- get_params.pkfit(obj)

   if(!("Rblood2plasma" %in% names(params))){
     params$Rblood2plasma <- 1
   }

    #call the model function with the params
    if(is.null(newdata)){
      newdata <- obj$data
    }

    preds <- tryCatch(do.call(modelfun,
            args = list(params = params,
                        time = newdata$Time,
                        dose = newdata$Dose,
                        iv.dose = newdata$iv,
                        medium = newdata$Media)),
            error = function(err) rep(NA_real_,
                                      nrow(newdata)))
    return(preds)

}

#' Extract fitted values from a `pkfit` object
#'
#' @param obj The pkfit object for which to extract the fitted values
#' @return A numeric vector of fitted values with length equal to
#'   `nrow(obj$data)`
#' @export
#' @author Caroline Ring
fitted.pkfit <- function(obj){
  predict.pkfit(obj)
}

#' Get residuals from a `pkfit` object
#'
#' @param obj A `pkfit` object
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), `Media` ("blood" or "plasma"), `Conc` (the
#'   observed concentration). Optionally, it may also include a variable
#'   `Detect`, which is TRUE if the observed concentration was above the limit
#'   of quantification and FALSE otherwise, and a variable `LOQ`, denoting the
#'   relevant limit of quantification. If there is no `Detect` variable, then
#'   all observed values are assumed to be detected. If an observed value is
#'   non-detect, then `Conc` should *not* contain NAs. Instead, it should
#'   already contain the desired substitution (e.g. `Conc` should be equal to
#'   LOQ, LOQ/2, LOQ/sqrt(2), 0, etc. -- whatever your desired substitution may
#'   be). But see argument `match_nondetect` for another option on how to
#'   compute residuals for non-detects. If concentration was reported as a mean
#'   and standard deviation of a group of subjects, then `Conc` should contain
#'   the reported mean, and `newdata` should also include variables `Conc_SD`
#'   (the reported standard deviation) and `N_Subjects` (the number of subjects
#'   in the group). If there are no variables `Conc_SD` and `N_Subjects`, then
#'   it is assumed that each observed concentration is for a single subject.
#' @param dose_norm Optional: TRUE to apply dose normalization before computing
#'   residuals (and before any log transformation); FALSE otherwise. Default
#'   NULL, to apply the transformation used for fitting
#'   (`obj$data_trans$fit_conc_dose`). If `TRUE`, then both observations and
#'   predictions will be divided by dose before any log-transformation is done.
#' @param log_trans Optional: TRUE to compute log-scale residuals (log(obs) -
#'   log(pred)); FALSE to compute natural-scale residuals (obs - pred). Default
#'   NULL, to apply the transformation used for fitting (
#'   `obj$data_trans$fit_log_conc`).
#' @param match_nondetect Logical. If TRUE (the default), then when an observed
#'   value was non-detect and the corresponding model-predicted value is below
#'   the reported LOQ, then the corresponding residual is reported as 0. If
#'   FALSE, then non-detect observed values will be treated like any other,
#'   using whatever substitution is already present in variable `Conc` of
#'   `newdata`.
#' @return Numeric: A vector of residuals for the specified `pkfit` object,
#'   data, and residual transformations.
#' @export
#' @author Caroline Ring

residuals.pkfit <- function(obj,
                            newdata = NULL,
                            dose_norm = NULL,
                            log_trans = NULL,
                            match_nondetect = TRUE){
  if(is.null(newdata)){
    newdata <- obj$data
  }
  #first get predictions
  preds <- predict.pkfit(obj = obj,
                   newdata = newdata)
  preds_orig <- preds
  #then get observations
if(is.null(dose_norm)){
  dose_norm <- obj$data_trans$fit_conc_dose
}

  if(dose_norm %in% FALSE){
    obs <- newdata$Conc
    obs_sd <- newdata$Conc_SD
  }else{
    obs <- newdata$Conc_Dose
    preds <- preds/newdata$Dose
    obs_sd <- newdata$Conc_SD/newdata$Dose
  }

  #log-transform if requested
  if(is.null(log_trans)){
    log_trans <- obj$data_trans$fit_log_conc
  }

  if("Conc_SD" %in% names(newdata)){
    if(log_trans %in% TRUE){
      tmplist <- convert_summary_to_log(sample_mean = obs,
                                        sample_SD = obs_sd)
      obs <- tmplist$logmean
      obs_sd <- tmplist$logSD
      preds <- log(preds)
    }
  }else{
if(log_trans %in% TRUE){
  obs <- log(obs)
  preds <- log(preds)
}
  }

  if(match_nondetect %in% TRUE){
    if(!("LOQ" %in% names(newdata))){
      newdata$LOQ <- 0
    }
    if(!("Detect" %in% names(newdata))){
      newdata$Detect <- TRUE
    }
    #set preds == obs if both are < LOQ
    preds <- ifelse(newdata$Detect %in% FALSE &
                      preds_orig < newdata$LOQ,
                    obs,
                    preds)
  }

  #calculate residuals
  resids <- obs - preds

  return(resids)

}

#' Alias to residuals.pkfit
#' @param obj A `pkfit` object
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), `Media` ("blood" or "plasma"), `Conc` (the
#'   observed concentration). Optionally, it may also include a variable
#'   `Detect`, which is TRUE if the observed concentration was above the limit
#'   of quantification and FALSE otherwise, and a variable `LOQ`, denoting the
#'   relevant limit of quantification. If there is no `Detect` variable, then
#'   all observed values are assumed to be detected. If an observed value is
#'   non-detect, then `Conc` should *not* contain NAs. Instead, it should
#'   already contain the desired substitution (e.g. `Conc` should be equal to
#'   LOQ, LOQ/2, LOQ/sqrt(2), 0, etc. -- whatever your desired substitution may
#'   be). But see argument `match_nondetect` for another option on how to
#'   compute residuals for non-detects. If concentration was reported as a mean
#'   and standard deviation of a group of subjects, then `Conc` should contain
#'   the reported mean, and `newdata` should also include variables `Conc_SD`
#'   (the reported standard deviation) and `N_Subjects` (the number of subjects
#'   in the group). If there are no variables `Conc_SD` and `N_Subjects`, then
#'   it is assumed that each observed concentration is for a single subject.
#' @param dose_norm Optional: TRUE to apply dose normalization before computing
#'   residuals (and before any log transformation); FALSE otherwise. Default
#'   NULL, to apply the transformation used for fitting
#'   (`obj$data_trans$fit_conc_dose`). If `TRUE`, then both observations and
#'   predictions will be divided by dose before proceeding.
#' @param log_trans Optional: TRUE to compute log-scale residuals (log(obs) -
#'   log(pred)); FALSE to compute natural-scale residuals (obs - pred). Default
#'   NULL, to apply the transformation used for fitting (
#'   `obj$data_trans$fit_log_conc`).
#' @param match_nondetect Logical. If TRUE (the default), then when an observed value was
#'   non-detect and the corresponding model-predicted value is below the
#'   reported LOQ, then the corresponding residual is reported as 0. If FALSE,
#'   then non-detect observed values will be treated like any other, using
#'   whatever substitution is already present in variable `Conc` of `newdata`.
#' @return Numeric: A vector of residuals for the specified `pkfit`
#'   object, data, and residual transformations.
#' @export
#' @author Caroline Ring
resids.pkfit <- function(...){
  residuals.pkfit(...)
}

#' Calculate root mean squared error (RMSE)
#'
#' @param obj A `pkfit` object
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), `Media` ("blood" or "plasma"), `Conc` (the
#'   observed concentration). Optionally, it may also include a variable
#'   `Detect`, which is TRUE if the observed concentration was above the limit
#'   of quantification and FALSE otherwise, and a variable `LOQ`, denoting the
#'   relevant limit of quantification. If there is no `Detect` variable, then
#'   all observed values are assumed to be detected. If an observed value is
#'   non-detect, then `Conc` should *not* contain NAs. Instead, it should
#'   already contain the desired substitution (e.g. `Conc` should be equal to
#'   LOQ, LOQ/2, LOQ/sqrt(2), 0, etc. -- whatever your desired substitution may
#'   be). But see argument `match_nondetect` for another option on how to
#'   compute residuals for non-detects. If concentration was reported as a mean
#'   and standard deviation of a group of subjects, then `Conc` should contain
#'   the reported mean, and `newdata` should also include variables `Conc_SD`
#'   (the reported standard deviation) and `N_Subjects` (the number of subjects
#'   in the group). If there are no variables `Conc_SD` and `N_Subjects`, then
#'   it is assumed that each observed concentration is for a single subject.
#' @param dose_norm Optional: TRUE to apply dose normalization before computing
#'   residuals (and before any log transformation); FALSE otherwise. Default
#'   NULL, to apply the transformation used for fitting
#'   (`obj$data_trans$fit_conc_dose`). If `TRUE`, then both observations and
#'   predictions will be divided by dose before proceeding.
#' @param log_trans Optional: TRUE to compute log-scale residuals (log(obs) -
#'   log(pred)); FALSE to compute natural-scale residuals (obs - pred). Default
#'   NULL, to apply the transformation used for fitting (
#'   `obj$data_trans$fit_log_conc`).
#' @param match_nondetect Logical. If TRUE (the default), then when an observed value was
#'   non-detect and the corresponding model-predicted value is below the
#'   reported LOQ, then the corresponding residual is reported as 0. If FALSE,
#'   then non-detect observed values will be treated like any other, using
#'   whatever substitution is already present in variable `Conc` of `newdata`.
#' @return Numeric scalar: the root mean squared error for the specified `pkfit`
#'   object, data, and residual transformations.
#' @export
#' @author Caroline Ring
rmse.pkfit <- function(obj,
                       newdata = NULL,
                       dose_norm = NULL,
                       log_trans = NULL,
                       match_nondetect = TRUE
){
  info_trans <- transform_obs_preds.pkfit(obj = obj,
                                          newdata = newdata,
                                          dose_norm = dose_norm,
                                          log_trans = log_trans,
                                          match_nondetect = match_nondetect)

   mse <- with(info_trans,
                (1/sum(obs_n)) * sum(
     (obs_n - 1) * obs_sd^2 + obs_n * obs^2 +
       -2 * preds * obs_n * obs  +
       obs_n * preds^2
   )
   )

   return(sqrt(mse))
}

#' Calculate R-squared for a fit
#'
#' @param obj A `pkfit` object
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), `Media` ("blood" or "plasma"), `Conc` (the
#'   observed concentration). Optionally, it may also include a variable
#'   `Detect`, which is TRUE if the observed concentration was above the limit
#'   of quantification and FALSE otherwise, and a variable `LOQ`, denoting the
#'   relevant limit of quantification. If there is no `Detect` variable, then
#'   all observed values are assumed to be detected. If an observed value is
#'   non-detect, then `Conc` should *not* contain NAs. Instead, it should
#'   already contain the desired substitution (e.g. `Conc` should be equal to
#'   LOQ, LOQ/2, LOQ/sqrt(2), 0, etc. -- whatever your desired substitution may
#'   be). But see argument `match_nondetect` for another option on how to
#'   compute residuals for non-detects. If concentration was reported as a mean
#'   and standard deviation of a group of subjects, then `Conc` should contain
#'   the reported mean, and `newdata` should also include variables `Conc_SD`
#'   (the reported standard deviation) and `N_Subjects` (the number of subjects
#'   in the group). If there are no variables `Conc_SD` and `N_Subjects`, then
#'   it is assumed that each observed concentration is for a single subject.
#' @param dose_norm Optional: TRUE to apply dose normalization before computing
#'   residuals (and before any log transformation); FALSE otherwise. Default
#'   NULL, to apply the transformation used for fitting
#'   (`obj$data_trans$fit_conc_dose`). If `TRUE`, then both observations and
#'   predictions will be divided by dose before proceeding.
#' @param log_trans Optional: TRUE to compute log-scale residuals (log(obs) -
#'   log(pred)); FALSE to compute natural-scale residuals (obs - pred). Default
#'   NULL, to apply the transformation used for fitting (
#'   `obj$data_trans$fit_log_conc`).
#' @param match_nondetect Logical. If TRUE (the default), then when an observed
#'   value was non-detect and the corresponding model-predicted value is below
#'   the reported LOQ, then the corresponding residual is reported as 0. If
#'   FALSE, then non-detect observed values will be treated like any other,
#'   using whatever substitution is already present in variable `Conc` of
#'   `newdata`.
#' @return Numeric scalar: the squared Pearson correlation coefficient between
#'   observed and predicted values for the specified `pkfit` object, data, and
#'   residual transformations.
#' @export
#' @author Caroline Ring
rsq.pkfit <- function(obj,
                      newdata = NULL,
                      dose_norm = NULL,
                      log_trans = NULL,
                      match_nondetect = TRUE){

info_trans <- transform_obs_preds.pkfit(obj = obj,
                          newdata = newdata,
                          dose_norm = dose_norm,
                          log_trans = log_trans,
                          match_nondetect = match_nondetect)
rsq <- with(info_trans,
     {
  #grand mean of predictions
  pred_bar <- sum(obs_n*preds)/sum(obs_n)
  #grand mean of observations
  x_bar <- sum(obs_n*obs)/sum(obs_n)
  #Calculate numerator of Pearson correlation coefficient
  r_num <- sum(preds * obs_n * obs) -
    pred_bar * sum(obs_n * obs) -
    x_bar * sum(obs_n * preds) +
    x_bar * pred_bar * sum(obs_n)
  #Calculate first denominator term of Pearson correlation coefficient
  r_denom_1 <- sqrt(sum(
    (obs_n - 1) * obs_sd^2 +
      obs_n * obs^2
  ) -
    2 * x_bar * sum(obs_n * obs) +
    sum(obs_n) * x_bar^2)
  #Calculate second denominator term of Pearson correlation coefficient
  r_denom_2 <- sqrt(sum(obs_n * preds^2) -
                      2 * pred_bar * sum(obs_n * preds) +
                      sum(obs_n) * pred_bar^2)

  #Put the pieces together to get Pearson correlation coefficient
  if((r_denom_1 > 0) %in% TRUE  &
     ( r_denom_2 > 0) %in% TRUE){
    r <- r_num /(r_denom_1 * r_denom_2)
  }else{
    r <- NA_real_
  }

  #Square it
   r^2
     })

  return(rsq)
}

#' Estimate the overall residual error standard deviation for a fit
#'
#' @param obj A `pkfit` object
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), `Media` ("blood" or "plasma"), `Conc` (the
#'   observed concentration). Optionally, it may also include a variable
#'   `Detect`, which is TRUE if the observed concentration was above the limit
#'   of quantification and FALSE otherwise, and a variable `LOQ`, denoting the
#'   relevant limit of quantification. If there is no `Detect` variable, then
#'   all observed values are assumed to be detected. If an observed value is
#'   non-detect, then `Conc` should *not* contain NAs. Instead, it should
#'   already contain the desired substitution (e.g. `Conc` should be equal to
#'   LOQ, LOQ/2, LOQ/sqrt(2), 0, etc. -- whatever your desired substitution may
#'   be). But see argument `match_nondetect` for another option on how to
#'   compute residuals for non-detects. If concentration was reported as a mean
#'   and standard deviation of a group of subjects, then `Conc` should contain
#'   the reported mean, and `newdata` should also include variables `Conc_SD`
#'   (the reported standard deviation) and `N_Subjects` (the number of subjects
#'   in the group). If there are no variables `Conc_SD` and `N_Subjects`, then
#'   it is assumed that each observed concentration is for a single subject.
#' @param dose_norm Optional: TRUE to apply dose normalization before computing
#'   residuals (and before any log transformation); FALSE otherwise. Default
#'   NULL, to apply the transformation used for fitting
#'   (`obj$data_trans$fit_conc_dose`). If `TRUE`, then both observations and
#'   predictions will be divided by dose before proceeding.
#' @param log_trans Optional: TRUE to compute log-scale residuals (log(obs) -
#'   log(pred)); FALSE to compute natural-scale residuals (obs - pred). Default
#'   NULL, to apply the transformation used for fitting (
#'   `obj$data_trans$fit_log_conc`).
#' @param optimx_args Optional: A named list of arguments to `optimx` (e.g.
#'   `method`, `itnmax`, and any `control` arguments). Default NULL to use the
#'   arguments specified in `obj$fit_args` and `obj$fit_control`.
#' @return A data.frame with three variables: `sigma`, the estimated overall
#'   error standard deviation for the specified `pkfit` object, data, and
#'   residual transformations; `convcode`, the convergence code returned by the
#'   optimizer; and `msg`, the error message returned by the optimizer if any.
#' @export
#' @author Caroline Ring
overall_error_sd.pkfit <- function(obj,
                                   newdata = NULL,
                                   dose_norm = NULL,
                                   log_trans = NULL,
                                   optimx_args = NULL,
                                   match_nondetect = TRUE){
  if(is.null(newdata)){
    newdata <- obj$data
  }

  if(is.null(log_trans)){
    log_trans <- obj$data_trans$fit_log_conc
  }

  if(is.null(dose_norm)){
    dose_norm <- obj$data_trans$fit_conc_dose
  }

  if(is.null(optimx_args)){
    optimx_args <- c(obj$fit_args,
                     list("control" = obj$fit_control))
  }

  resids <- resids.pkfit(obj = obj,
                         newdata = newdata,
                         log_trans = log_trans,
                         dose_norm = dose_norm,
                         match_nondetect = match_nondetect)

  sd_start <- sd(resids)

  params <- get_params.pkfit(obj)
  if(!("Rblood2plasma" %in% names(params))){
    params$Rblood2plasma <- 1
  }
  model_params <- get_model_paramnames(model = obj$model_info$model)
  params <- params[names(params) %in% model_params]

  #Hold model params constant & optimize only a pooled sigma
  sigma_fit <- tryCatch(do.call(optimx::optimx,
                                args = c(
                                  #
                                  list(par = c("sigma" = sd_start),
                                       fn = log_likelihood,
                                       lower = c("sigma" = 1e-4),
                                       upper = c("sigma" = 1e4)),
                                  #method and control
                                  optimx_args,
                                  #other args to log_likelihood()
                                  list(
                                    const_params = params, #fitted model params
                                    DF = newdata,
                                    modelfun = obj$model_info$model.type,
                                    model = obj$model_info$model,
                                    fit_conc_dose = dose_norm,
                                    fit_log_conc = log_trans,
                                    force_finite = (optimx_args$method %in% "L-BFGS-B"),
                                    negative = TRUE
                                  )
                                )
  ),
  error = function(err) {
    return(data.frame(sigma = NA_real_,
                      convcode = -99,
                      msg = err$message))
  })

  if(!("msg" %in% names(sigma_fit))){
    sigma_fit$msg <- "none"
  }

  return(sigma_fit[c("sigma",
                      "convcode",
                      "msg")])

}

#' Calculate summary TK statistics from a fitted object
#'
#' @param obj A `pkfit` object
calc_TK_stats.pkfit <- function(obj){
  model <- obj$model_info$model

  expected_names <- c(get_model_paramnames(model),
                      paste(get_model_paramnames(model),
                            "sd",
                            sep = "_")
  )
  DF <- data.frame(as.list(rep(NA_real_,
                               length(expected_names))))
  names(DF) <- expected_names
  params <- get_params.pkfit(obj)
  params_sd <- get_params_sd.pkfit(obj)

  for(this_param in intersect(names(DF),
                               names(params))){
    DF[this_param] <- params[this_param]
  }

  for(this_param_sd in intersect(names(DF),
                              names(params_sd))){
    DF[this_param_sd] <- params_sd[this_param_sd]
  }

  if(model %in% "flat"){
    expected_names <- c(get_model_paramnames(model),
                        paste(get_model_paramnames(model),
                              "sd",
                              sep = "_")
    )

    #initalize by filling with NA_real_
    DF[setdiff(expected_names,
               names(DF))] <- NA_real_
    #sort all the columns in a standard expected order
    DF <- DF[expected_names]
  }else if(model %in% "1compartment"){
    expected_names <- c(get_model_paramnames(model),
                        paste(get_model_paramnames(model),
                              "sd",
                              sep = "_"),
                        "halflife",
                        "CLtot",
                        "Css_iv_1mgkg",
                        "AUC_inf_iv_1mgkg",
                        "CLtot_Fgutabs",
                        "Css_oral_1mgkg",
                        "tmax_oral",
                        "Cmax_oral_1mgkg",
                        "AUC_inf_oral_1mgkg")
#initalize by filling with NA_real_
    DF[setdiff(expected_names,
               names(DF))] <- NA_real_
    #sort all the columns in a standard expected order
    DF <- DF[expected_names]

    #if Fgutabs and Vdist were computed separately, create a temp Fgutabs_Vdist
    if(!is.na(DF$Fgutabs)){
    DF$Fgutabs_Vdist <- DF$Fgutabs/DF$Vdist
    }

    #half-life, tmax, Cmax:
    #see https://www.boomer.org/c/p4/c08/c0803.php

    #half-life
    DF$halflife <-log(2) / DF$kelim

    if(!(is.na(DF$Vdist))){
      #CLtot:
      #If kelim and Vdist are available, CLtot = kelim * Vdist
      DF$CLtot <- DF$kelim * DF$Vdist

      #Css for 1 mg/kg/day IV infusion dose (if Vdist available)
      DF$Css_iv_1mgkg <- 1/(24 * DF$CLtot)

      #AUC_infinity for IV bolus dose of 1 mg/kg
      DF$AUC_inf_iv_1mgkg <- auc_1comp(params = list(
        "kelim" = DF$kelim,
        "Vdist" = DF$Vdist
      ),
      time = Inf,
      dose = 1,
      iv.dose = TRUE,
      medium = "plasma")
    }

    if(!is.na(DF$Fgutabs_Vdist)){
    #if kelim and Fgutabs/Vdist are available, can only get CLtot/Fgutabs
    #kelim * (Vdist/Fgutabs) = CLtot/Fgutabs
    DF$CLtot_Fgutabs <- DF$kelim / DF$Fgutabs_Vdist

    #Css_for a 1 mg/kg/day oral infusion dose
    DF$Css_oral_1mgkg <- DF$Fgutabs_Vdist * (1/(24*DF$kelim))

    #tmax -- only available if kgutabs was fitted
    DF$tmax_oral <-log(DF$kgutabs / DF$kelim) /
      (DF$kgutabs - DF$kelim)

    #Cmax for oral bolus dose of 1 mg/kg
    DF$Cmax_oral_1mgkg <- cp_1comp(params = list(
      "kelim" = DF$kelim,
      "Fgutabs_Vdist" = DF$Fgutabs_Vdist,
      "kgutabs" = DF$kgutabs
    ),
    time = DF$tmax_oral,
    dose = 1,
    iv.dose = FALSE,
    medium = "plasma")

    #AUC_infinity for oral bolus dose of 1 mg/kg
    DF$AUC_inf_oral_1mgkg <- auc_1comp(params = list(
      "kelim" = DF$kelim,
      "Fgutabs_Vdist" = DF$Fgutabs_Vdist,
      "kgutabs" = DF$kgutabs
    ),
    time = Inf,
    dose = 1,
    iv.dose = FALSE,
    medium = "plasma")
    }

    #switch back to original Fgutabs_Vdist (i.e., NA if Fgutabs and Vdist were
    #fitted separately)
   if(!is.na(DF$Fgutabs)){
     DF$Fgutabs_Vdist <- NA_real_
   }
    #ensure that variables are in a standard expected order
    DF <- DF[expected_names]
  }else if(model %in% "2compartment"){
    expected_names <- c(get_model_paramnames(model),
                        paste(get_model_paramnames(model), "sd", sep = "_"),
                        "alphabeta_sum",
                        "alphabeta_prod",
                        "alpha",
                        "beta",
                        "halflife_beta",
                        "halflife_alpha",
                        "CLtot",
                        "A_iv_1mgkg",
                        "B_iv_1mgkg",
                        "Vss",
                        "Vbeta",
                        "Css_iv_1mgkg",
                        "AUC_inf_iv_1mgkg",
                        "CLtot_Fgutabs",
                        "Vss_Fgutabs",
                        "A_oral_1mgkg",
                        "B_oral_1mgkg",
                        "Vbeta_Fgutabs",
                        "Css_oral_1mgkg",
                        "halflife_abs",
                        "tmax_oral",
                        "Cmax_oral_1mgkg",
                        "AUC_inf_oral_1mgkg")

    #initalize by filling with NA_real_
    DF[setdiff(expected_names,
               names(DF))] <- NA_real_
    #sort all the columns in a standard expected order
    DF <- DF[expected_names]

    #in case Fgutabs and V1 were fitted separately, compute Fgutabs/V1
    if(!(is.na(DF$Fgutabs))){
    DF$Fgutabs_V1 <-  DF$Fgutabs/DF$V1
    }

    #compute A, B, alpha, beta
    #see https://www.boomer.org/c/p4/c19/c1902.php
    DF$alphabeta_sum <- DF$kelim + DF$k12 + DF$k21
    DF$alphabeta_prod <- DF$kelim * DF$k21
    DF$alpha <- (DF$alphabeta_sum + sqrt(DF$alphabeta_sum^2 - 4*DF$alphabeta_prod))/2
    DF$beta <- (DF$alphabeta_sum - sqrt(DF$alphabeta_sum^2 - 4*DF$alphabeta_prod))/2
    #half-lives for each phase
    DF$halflife_beta <- log(2) / DF$beta
    DF$halflife_alpha <- log(2) / DF$alpha

    #if there was IV data and V1 was fitted:
    if(!(is.na(DF$V1))){
    #Total clearance
    DF$CLtot <- DF$kelim * DF$V1
    #A and B for IV data (for dose of 1 mg/kg)
    DF$A_iv_1mgkg <- 1*(DF$alpha - DF$k21) / (DF$V1 * (DF$alpha - DF$beta))
    DF$B_iv_1mgkg<- 1*(DF$k21 - DF$beta) / (DF$V1 * (DF$alpha - DF$beta))

    #Vbeta
    #Terminal volume of distribution
    #see https://www.boomer.org/c/p4/c19/c1905.php
    DF$Vbeta <- DF$V1 * DF$kelim / DF$beta

    #Vss
    #apparent volume of distribution at steady state
    #see https://www.boomer.org/c/p4/c19/c1905.php
    DF$Vss <- DF$V1 * (DF$k21 + DF$k12) / DF$k21

    #Get Css = average plasma concentration for 1 mg/kg/day every 1 days
    DF$Css_iv_1mgkg <- 1/(24 * DF$CLtot)

    #AUC at t = infinity
    DF$AUC_inf_iv_1mgkg <- auc_2comp(params = list("kelim" = DF$kelim,
                                                   "V1" = DF$V1,
                                                   "k12" = DF$k12,
                                                   "k21" = DF$k21),
                                     time = Inf,
                                     dose = 1,
                                     iv.dose = TRUE,
                                     medium = "plasma")

    }

    if(!(is.na(DF$Fgutabs_V1))){
      #Total clearance divided by Fgutabs, in case only Fgutabs/V1 was available
      DF$CLtot_Fgutabs <- DF$kelim / DF$Fgutabs_V1
      #Vss/Fgutabs, in case only Fgutabs/V1 was available and not V1 by itself
      DF$Vss_Fgutabs <- (1/DF$Fgutabs_V1) * (DF$k21 + DF$k12) / DF$k21

      #A and B when oral data were available (for dose of 1 mg/kg)
      DF$A_oral_1mgkg <- (DF$kgutabs * DF$Fgutabs_V1 *
                       (DF$alpha - DF$k21)) /
        ( (DF$kgutabs - DF$alpha) * (DF$alpha - DF$beta))
      DF$B_oral_1mgkg <- (DF$kgutabs * DF$Fgutabs_V1 *
                       (DF$k21 - DF$beta)) /
        ( (DF$kgutabs - DF$beta) * (DF$alpha - DF$beta))
      DF$Vbeta_Fgutabs <- (1/DF$Fgutabs_V1) * DF$kelim / DF$beta
      DF$Css_oral_1mgkg <- DF$Fgutabs_V1 * (1/(24*DF$kelim))
      DF$halflife_abs <- log(2) / DF$kgutabs

    #tmax for oral dose
    #I don't think it can be done analytically
    #so do it numerically
    #use uniroot() to search for zeros of time deriv of Cp vs. time
    #for each set of parameters
    #look between 0 and what would be the 1-compartment tmax
    #extending the interval if necessary
    #if uniroot() fails, then just return NA
DF$tmax_oral <- tryCatch(
      uniroot( f = function(x){
        cp_2comp_dt(params = list("kelim" = DF$kelim,
                                  "Fgutabs_V1" = DF$Fgutabs_V1,
                                  "kgutabs" = DF$kgutabs,
                                  "k12" = DF$k12,
                                  "k21" = DF$k21),
                    time = x,
                    dose = 1,
                    iv.dose = FALSE,
                    medium = "plasma")
      },
      lower = 0,
      upper = log(DF$kgutabs / DF$kelim) / (DF$kgutabs - DF$kelim),
      extendInt = "downX", #function should be decreasing
      maxiter = 1000,
      tol = .Machine$double.eps)$root,
      error = function(err) return(NA_real_))

    #Cmax for oral bolus dose of 1 mg/kg
DF$Cmax_oral_1mgkg <- cp_2comp(params = list("kelim" = DF$kelim,
                                                       "Fgutabs_V1" = DF$Fgutabs_V1,
                                                       "kgutabs" = DF$kgutabs,
                                                       "k12" = DF$k12,
                                                       "k21" = DF$k21),
                                         time = DF$tmax_oral,
                                         dose = 1,
                                         iv.dose = FALSE,
                                         medium = "plasma")


    #AUC_infinity for oral bolus dose of 1 mg/kg
DF$AUC_inf_oral_1mgkg <- auc_2comp(params = list("kelim" = DF$kelim,
                                                           "Fgutabs_V1" = DF$Fgutabs_V1,
                                                           "kgutabs" = DF$kgutabs,
                                                           "k12" = DF$k12,
                                                           "k21" = DF$k21),
                                             time = Inf,
                                             dose = 1,
                                             iv.dose = FALSE,
                                             medium = "plasma")
    }else{
      DF$CLtot_Fgutabs <- NA_real_
      #Vss/Fgutabs, in case only Fgutabs/V1 was available and not V1 by itself
      DF$Vss_Fgutabs <- NA_real_

      #A and B when oral data were available (for dose of 1 mg/kg)
      DF$A_oral_1mgkg <- NA_real_
      DF$B_oral_1mgkg <- NA_real_
      DF$Vbeta_Fgutabs <- NA_real_
      DF$Css_oral_1mgkg <- NA_real_
      DF$halflife_abs <- NA_real_
      DF$tmax_oral <- NA_real_
      DF$Cmax_oral_1mgkg <- NA_real_
      DF$AUC_inf_oral_1mgkg <- NA_real_
}

if(!(is.na(DF$Fgutabs))){
  DF$Fgutabs_V1 <- NA_real_
}
    #ensure that variables are in a standard expected order
    DF <- DF[expected_names]
  } #end if model %in% "2compartment"


return(DF)
}

#' Compute fractional residuals
#'
#' Compute residuals as (pred/obs)
#' @param obj A `pkfit` object
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), `Media` ("blood" or "plasma"), `Conc` (the
#'   observed concentration). Optionally, it may also include a variable
#'   `Detect`, which is TRUE if the observed concentration was above the limit
#'   of quantification and FALSE otherwise, and a variable `LOQ`, denoting the
#'   relevant limit of quantification. If there is no `Detect` variable, then
#'   all observed values are assumed to be detected. If an observed value is
#'   non-detect, then `Conc` should *not* contain NAs. Instead, it should
#'   already contain the desired substitution (e.g. `Conc` should be equal to
#'   LOQ, LOQ/2, LOQ/sqrt(2), 0, etc. -- whatever your desired substitution may
#'   be). But see argument `match_nondetect` for another option on how to
#'   compute residuals for non-detects. If concentration was reported as a mean
#'   and standard deviation of a group of subjects, then `Conc` should contain
#'   the reported mean, and `newdata` should also include variables `Conc_SD`
#'   (the reported standard deviation) and `N_Subjects` (the number of subjects
#'   in the group). If there are no variables `Conc_SD` and `N_Subjects`, then
#'   it is assumed that each observed concentration is for a single subject.
#' @param dose_norm Optional: TRUE to apply dose normalization before computing
#'   residuals (and before any log transformation); FALSE otherwise. Default
#'   NULL, to apply the transformation used for fitting
#'   (`obj$data_trans$fit_conc_dose`). If `TRUE`, then both observations and
#'   predictions will be divided by dose before proceeding.
#' @param log_trans Optional: TRUE to compute log-scale residuals (log(obs) -
#'   log(pred)); FALSE to compute natural-scale residuals (obs - pred). Default
#'   NULL, to apply the transformation used for fitting (
#'   `obj$data_trans$fit_log_conc`).
#' @return A numeric vector of fractional residuals.
#' @export
#' @author Caroline Ring
frac_resids.pkfit <- function(obj,
                              newdata = NULL,
                              dose_norm = NULL,
                              log_trans = NULL,
                              match_nondetect = TRUE){

  info_trans <- transform_obs_preds.pkfit(obj = obj,
                                          newdata = newdata,
                                          dose_norm = dose_norm,
                                          log_trans = log_trans,
                                          match_nondetect = match_nondetect)

  frac_resids <- with(info_trans,
                      preds/obs)

  return(frac_resids)
}

#' Helper function to transform observations and predictions
#'
#' @param obj A `pkfit` object
#' @param newdata Optional: A `data.frame` of new data for which to make
#'   predictions. Default NULL, to use the fitted data (`obj$data`). If
#'   `newdata` is provided, it must be a `data.frame` containing variables named
#'   `Time` (in hours), `Dose` (in mg/kg), `iv` (TRUE for IV administration,
#'   FALSE for oral administration), `Media` ("blood" or "plasma"), `Conc` (the
#'   observed concentration). Optionally, it may also include a variable
#'   `Detect`, which is TRUE if the observed concentration was above the limit
#'   of quantification and FALSE otherwise, and a variable `LOQ`, denoting the
#'   relevant limit of quantification. If there is no `Detect` variable, then
#'   all observed values are assumed to be detected. If an observed value is
#'   non-detect, then `Conc` should *not* contain NAs. Instead, it should
#'   already contain the desired substitution (e.g. `Conc` should be equal to
#'   LOQ, LOQ/2, LOQ/sqrt(2), 0, etc. -- whatever your desired substitution may
#'   be). But see argument `match_nondetect` for another option on how to
#'   compute residuals for non-detects. If concentration was reported as a mean
#'   and standard deviation of a group of subjects, then `Conc` should contain
#'   the reported mean, and `newdata` should also include variables `Conc_SD`
#'   (the reported standard deviation) and `N_Subjects` (the number of subjects
#'   in the group). If there are no variables `Conc_SD` and `N_Subjects`, then
#'   it is assumed that each observed concentration is for a single subject.
#' @param dose_norm Optional: TRUE to apply dose normalization before computing
#'   residuals (and before any log transformation); FALSE otherwise. Default
#'   NULL, to apply the transformation used for fitting
#'   (`obj$data_trans$fit_conc_dose`). If `TRUE`, then both observations and
#'   predictions will be divided by dose before proceeding.
#' @param log_trans Optional: TRUE to compute log-scale residuals (log(obs) -
#'   log(pred)); FALSE to compute natural-scale residuals (obs - pred). Default
#'   NULL, to apply the transformation used for fitting (
#'   `obj$data_trans$fit_log_conc`).
#' @return A list with named components: "obs", containing the observations,
#'   "obs_sd", containing the observed standard deviations (filled with zeros if
#'   `newdata` does not have a variable `Conc_SD`), "obs_n" containing the
#'   number of subjects summarized in each observation (filled with 1
#'   if`newdata` does not have a variable `N_Subjects`) and "preds" containing
#'   the predictions. Observations, observed SDs, and predictions are
#'   transformed according to `dose_norm` and `log_trans`.
#' @author Caroline Ring
transform_obs_preds.pkfit <- function(obj,
                               newdata = NULL,
                               dose_norm = NULL,
                               log_trans = NULL,
                               match_nondetect = TRUE){
  if(is.null(newdata)){
    newdata <- obj$data
  }
  #first get predictions
  preds <- predict.pkfit(obj = obj,
                         newdata = newdata)
  preds_orig <- preds
  #then get observations
  if(is.null(dose_norm)){
    dose_norm <- obj$data_trans$fit_conc_dose
  }

  if(dose_norm %in% FALSE){
    obs <- newdata$Conc
    obs_sd <- newdata$Conc_SD
  }else{
    obs <- newdata$Conc_Dose
    preds <- preds/newdata$Dose
    obs_sd <- newdata$Conc_SD/newdata$Dose
  }

  #log-transform if requested
  if(is.null(log_trans)){
    log_trans <- obj$data_trans$fit_log_conc
  }

  if("Conc_SD" %in% names(newdata)){
    if(log_trans %in% TRUE){
      tmplist <- convert_summary_to_log(sample_mean = obs,
                                        sample_SD = obs_sd)
      obs <- tmplist$logmean
      obs_sd <- tmplist$logSD
      preds <- log(preds)
    }
  }else{
    obs_sd <- 0
    if(log_trans %in% TRUE){
      obs <- log(obs)
      preds <- log(preds)
    }
  }

  if("N_Subjects" %in% names(newdata)){
    obs_n <- newdata$N_Subjects
  }else{
    obs_n <- 1
  }

  if(match_nondetect %in% TRUE){
    #set preds == obs
    preds <- ifelse(newdata$Detect %in% FALSE &
                      preds_orig < newdata$LOQ,
                    obs,
                    preds)
  }

  return(list("obs" = obs,
              "obs_sd" = obs_sd,
              "obs_n" = obs_n,
              "preds" = preds))
}



#' Plot data and fit for a `pkfit` object.
#'
#' @param obj The `pkfit` object containing the fit to be plotted.
#' @param newdata Optional: A new concentration vs. time data set to plot,
#'   instead of `obj$data`.
#' @param dose_norm Optional: TRUE to plot dose-normalized concentrations; FALSE
#'   to plot non-dose normalized concentrations. Default NULL, to use the
#'   transformation applied during the fit (`obj$data_trans$fit_conc_dose`.)
#' @param plot_log10_conc Optional: TRUE to apply log10-scaling to the y axis
#'   (concentration or dose-normalized concentration). FALSE to leave the y axis
#'   on the natural scale. Default NULL, to use the transformation applied
#'   during the fit (`obj$data_trans$fit_log_conc`.)
#' @param n_interp_time Number of time points to interpolate between each unique
#'   time point in the data, in order to plot a smoother curve for the model
#'   fit. Default 10.
#' @return A `ggplot2` plot object.
plot.pkfit <- function(obj,
                       newdata = NULL,
                       plot_dose_norm = NULL,
                       plot_log10_conc = NULL,
                       n_interp_time = 10){
 if(is.null(newdata)){
   newdata <- obj$data
 }

  if(is.null(dose_norm)){
    plot_dose_norm <- obj$data_trans$fit_conc_dose
  }

  if(is.null(plot_log10_conc)){
    plot_log10_conc <- obj$data_trans$fit_log_conc
  }


pred_DF <- get_predDF.pkfit(obj = obj,
                            newdata = newdata,
                            dose_norm = plot_dose_norm,
                            n_interp_time = n_interp_time)

 #generate plot title
 plot_title <- paste0(obj$data_info$DTXSID,
                      " (", unique(obj$data$Compound), ")\n",
                      "Species = ", obj$data_info$Species, ", ",
                      "Doses = ", paste(signif(
                        sort(unique(obj$data$Dose)),
                        3),
                        collapse = ", "), " mg/kg\n",
                      "Analysis Type = ", obj$analysis_type, "\n",
                      "Fitting options: ",
                      "log-transform ",
                      obj$data_trans$fit_log_conc,
                      "; dose-normalize ",
                      obj$data_trans$fit_conc_dose,
                      "; rescale time ",
                      obj$data_trans$rescale_time)

 #generate plot subtitles

 #get coeffs of model
 par <- coef.pkfit(obj)
 #paste into a comma-separated list
 par_char <- paste(
   paste(names(par),
         signif(par, 3), #keep 3 sigfigs
         sep=" = "),
   collapse = ", ")

 plot_subtitle <- paste0("Model = ",
                        obj$model_info$model)
 plot_subtitle <- paste0(plot_subtitle, "\n",
                         par_char)

 #plot concentration vs. time data

 #note that any pkfit object is also a pkdata
 #object, so we can just call plot.pkdata()
 p <- plot.pkdata(obj = obj,
                      newdata = newdata,
                      plot_dose_norm = plot_dose_norm,
                      plot_log10_conc = plot_log10_conc)
#add the model-fit curve(s)
 p <- p + geom_line(pred_DF) +
   #replace the existing title with the new title/subtitle
   labs(title = plot_title,
        subtitle = plot_subtitle) +
 return(p)
}

#' Produce a summary data.frame for a `pkfit` fitted object.
#'
#' @param obj The `pkfit` object to summarize.
#' @param dose_norm Optional: TRUE to calculate metrics using dose-normalized
#'   residuals, FALSE to calculate metrics using non-dose-normalized residuals.
#'   Default NULL to use the transformation applied during the fit
#'   (`obj$data_trans$fit_conc_dose`).
#' @param log_trans Optional: TRUE to calculate metrics using log-scale
#'   residuals, FALSE to calculate metrics using natural-scale residuals.
#'   Default NULL to use the transformation applied during the fit
#'   (`obj$data_trans$fit_log_conc`).
#' @return A named list with the following items:
#' * "data_info": A `data.frame` of `obj$data_info`
#' * "data_trans": A `data.frame` of `obj$data_trans`
#' * "resid_trans": A `data.frame` of `dose_norm` and `log_trans`, in case they are different from `obj$data_trans`
#' * "gof": A `data.frame` of goodness-of-fit metrics with one row and three variables: `AIC` (the Akaike Information Criterion), `BIC` (the Bayesian Information Criterion), and `LogLikelihood`.
#' * "resids": A `data.frame` of metrics calculated on residuals, including RMSE, r-squared, overall error standard deviation, and summary statistics on fractional residuals (mean, min, max, and quantiles).
#' * "params": A `data.frame` of fitted model parameters, with one row and as many variables as there are model parameters.
#' * "params_sd": A `data.frame` of fitted model parameter standard deviations, with one row and as many variables as there are model parameters.
#' * "tk_stats": A `data.frame` of predicted toxicokinetic statistics as from [calc_TK_stats.pkfit()] (excluding the fitted parameters and SDs). Two additional quantities are added, `AUC_tlast_oral_1mgkg` and `AUC_tlast_iv_1mgkg`, which predict AUC at the latest time point in the data for oral administration and IV administration; they will be NA if `obj$data` contains no oral data or no IV data (respectively). respectively).
#' * "nca_stats": A `data.frame` of non-compartmental estimates of toxicokinetic statistics, as from [nca.pkfit()].
#' @export
#' @author Caroline Ring
summary.pkfit <- function(obj,
                          dose_norm = NULL,
                          log_trans = NULL,
                          match_nondetect = TRUE){
  if(is.null(dose_norm)){
    dose_norm <- obj$data_trans$fit_conc_dose
  }

  if(is.null(log_trans)){
    log_trans <- obj$data_trans$fit_log_conc
  }

 rmse <- rmse.pkfit(obj = obj,
                    dose_norm =  dose_norm,
                    log_trans = log_trans)

 rsq <- rsq.pkfit(obj = obj,
                  dose_norm =  dose_norm,
                  log_trans = log_trans,
                  match_nondetect = match_nondetect)
 overall_err_sd <- overall_error_sd.pkfit(obj = obj,
                                          dose_norm =  dose_norm,
                                          log_trans = log_trans,
                                          match_nondetect = match_nondetect)
 frac_resid <- frac_resids.pkfit(obj = obj,
                                dose_norm =  dose_norm,
                                log_trans = log_trans,
                                match_nondetect = match_nondetect)

 tk_stats <- calc_TK_stats.pkfit(obj = obj)

if(any(obj$data$Route %in% "po")){
  newdata <- subset(obj$data, Route %in% "po")
  AUC_tlast_oral_1mgkg <- predict.pkfit(obj = obj,
                                        newdata = data.frame(Time = max(newdata$Time),
                                                             Dose = 1,
                                                             iv = FALSE,
                                                             Media = "plasma"))
}else{
  AUC_tlast_oral_1mgkg <- NA_real_
}

 if(any(obj$data$Route %in% "iv")){
   newdata <- subset(obj$data, Route %in% "iv")
   AUC_tlast_iv_1mgkg <- predict.pkfit(obj = obj,
                                         newdata = data.frame(Time = max(newdata$Time),
                                                              Dose = 1,
                                                              iv = TRUE,
                                                              Media = "plasma"))
 }else{
   AUC_tlast_iv_1mgkg <- NA_real_
 }

 tk_stats$AUC_tlast_oral_1mgkg <- AUC_tlast_oral_1mgkg
 tk_stats$AUC_tlast_iv_1mgkg <- AUC_tlast_iv_1mgkg

 params <- as.data.frame(get_params.pkfit(obj))
 params_sd <- as.data.frame(get_params_sd.pkfit(obj))
 #remove parameters from tk_stats
 tk_stats <- tk_stats[setdiff(names(tk_stats),
                              c(names(params),
                                names(params_sd)
                                )
                              )
                      ]

 nca_stats <- nca.pkfit(obj = obj,
                        dose_norm = TRUE)

#Return a list of data.frames

 frac_resid_q <- as.list(quantile(frac_resid,
                          probs = c(0.01,0.05,0.1,
                                    0.25,0.5,0.75,
                                    0.9,0.95,0.99)))
 names(frac_resid_q) <- paste0("frac_resid.P",
                              100*c(0.01,0.05,0.1,
                                    0.25,0.5,0.75,
                                    0.9,0.95,0.99))

 out_list <- list(
   "data_info" = as.data.frame(obj$data_info),
   "data_trans" = as.data.frame(obj$data_trans),
   "resid_trans" = data.frame("dose_norm" = dose_norm,
                              "log_trans" = log_trans),
   "gof" = data.frame(
     AIC = AIC.pkfit(obj = obj,
                     newdata = NULL,
                     dose_norm = dose_norm,
                     log_trans = log_trans),
     BIC = BIC.pkfit(obj = obj,
                     newdata = NULL,
                     dose_norm = dose_norm,
                     log_trans = log_trans),
     LogLikelihood = loglike.pkfit(obj = obj,
                                   newdata = NULL,
                                   dose_norm = dose_norm,
                                   log_trans = log_trans)
   ),
   "resids" = as.data.frame(c(
     list(
       "RMSE" = rmse,
       "rsq" = rsq,
       "frac_resid.mean" = mean(frac_resid),
       "frac_resid.median" = median(frac_resid),
       "frac_resid.min" = min(frac_resid),
       "frac_resid.max" = max(frac_resid)
     ),
     frac_resid.q)),
   "params" = params,
   "params_sd" = params_sd,
   "tkstats" = tk_stats,
   "nca" = nca_stats
 )

 return(out_list)
}

AIC.pkfit <- function(obj,
                      newdata = NULL,
                      dose_norm = NULL,
                      log_trans = NULL){
#if no new data and no new data transformations,
  #just grab the fitted AIC
  if(is.null(newdata) &
     (is.null(dose_norm) &
      is.null(log_trans)) |
     (dose_norm %in% obj$data_trans$fit_conc_dose &
      log_trans %in% obj$data_trans$fit_log_conc)){
  AIC <- obj$fit_gof[["AIC"]]
  }else{
    #if new data and/or new data transformations, recalculate AIC
    if(is.null(newdata)){
      newdata <- obj$data
    }
    if(is.null(dose_norm)){
      dose_norm <- obj$data_trans$fit_conc_dose
    }
    if(is.null(log_trans)){
      log_trans <- obj$data_trans$fit_log_conc
    }
    params <- get_params.pkfit(obj)
    #evaluate log-likelihood
    loglike <- log_likelihood(params = params,
                              const_params = NULL,
                              DF = newdata,
                              fit_conc_dose = dose_norm,
                              fit_log_conc = log_trans,
                              force_finite = FALSE,
                              negative = FALSE)
    #calculate AIC
    AIC <- 2 * length(params) - 2 * loglike
  }
return(AIC)
}

loglike.pkfit <- function(obj,
                          newdata = NULL,
                          dose_norm = NULL,
                          log_trans = NULL){
  #if no new data and no new data transformations,
  #just grab the fitted log-likelihood
  if(is.null(newdata) &
     (is.null(dose_norm) &
      is.null(log_trans)) |
     (dose_norm %in% obj$data_trans$fit_conc_dose &
      log_trans %in% obj$data_trans$fit_log_conc)){
   loglike <- obj$fit_gof[["LogLikelihood"]]
  }else{
    #if new data and/or new data transformations, recalculate AIC
    if(is.null(newdata)){
      newdata <- obj$data
    }
    if(is.null(dose_norm)){
      dose_norm <- obj$data_trans$fit_conc_dose
    }
    if(is.null(log_trans)){
      log_trans <- obj$data_trans$fit_log_conc
    }
    params <- get_params.pkfit(obj)
    #evaluate log-likelihood
    loglike <- log_likelihood(params = params,
                              const_params = NULL,
                              DF = newdata,
                              fit_conc_dose = dose_norm,
                              fit_log_conc = log_trans,
                              force_finite = FALSE,
                              negative = FALSE)
  }
return(loglike)
}

#' Helper function for plotting curve fits to data
#'
#' Produce a data.frame of curve fits
#'
#' Create a data.frame that interpolates between time points of concentration
#' vs. time data, and gets model-predicted concentration for each combination of
#' route, medium, and (if `dose_norm == FALSE`), dose.
#'
#' @param obj
#' @return  A `data.frame` with variables `Time`, `Dose`, `Route`, `Medium`, and `Conc`,
get_predDF.pkfit <- function(obj,
                             newdata = NULL,
                             dose_norm = NULL,
                             n_interp_time = 10){
  if(is.null(newdata)){
    newdata <- obj$data
  }

  if(is.null(dose_norm)){
    dose_norm <- obj$data_trans$fit_conc_dose
  }

  #produce interpolated time points: n_interp_time points between every existing
  #time point
  timepoints <- sort(unique(c(0,newdata$Time)))
  time_interp <- sapply(1:(length(timepoints)-1),
                        function(i){
                          seq(from = timepoints[i],
                              to = timepoints[i+1],
                              length.out = n_interp_time)
                        }
  )

  #now repeat the time points for each unique combination of Study, Route,
  #Media, Dose in the data
  if(dose_norm %in% FALSE){
    #if not plotting dose-normalized,
    #make predictions for each dose in newdata
    combos <- unique(newdata[c("Route",
                               "iv",
                               "Media",
                               "Dose")])
  }else{
    #if plotting dose-normalized:
    #just make predictions for a dose of 1 mg/kg/day
    combos <- unique(newdata[c("Route",
                               "iv",
                               "Media")])
    combos$Dose <- 1
  }

  pred_DF <- combos[rep(seq_len(nrow(combos)),
                        each = length(timepoints))]
  pred_DF$Time <- rep(timepoints, nrow(combos))

  #now make the predictions
  pred_DF$Conc <- predict.pkfit(obj, newdata = pred_DF)

  #get predicted conc normalized by dose
  pred_DF$Conc_Dose<- pred_DF$Conc/pred_DF$Dose


  #create generic named columns for plotting, containing either dose-normalized or
  #non-dose-normalized concentrations, depending on `dose_norm`
  if(dose_norm %in% TRUE){
    pred_DF$conc_plot <- pred_DF$Conc/pred_DF$Dose
  }else{
    pred_DF$conc_plot <- pred_DF$Conc
  }

  return(pred_DF)
}

BIC.pkfit <- function(obj,
                      newdata = NULL,
                      dose_norm = NULL,
                      log_trans = NULL){
  #if no new data and no new data transformations,
  #just grab the fitted BIC
  if(is.null(newdata) &
     (is.null(dose_norm) &
      is.null(log_trans)) |
     (dose_norm %in% obj$data_trans$fit_conc_dose &
      log_trans %in% obj$data_trans$fit_log_conc)){
    k <- length(coef.pkfit(obj))
    n <- nrow(obj$data)
    BIC <- obj$fit_gof[["AIC"]] + k * log(n) - 2*k
  }else{
    #if new data and/or new data transformations, recalculate BIC
    if(is.null(newdata)){
      newdata <- obj$data
    }
    if(is.null(dose_norm)){
      dose_norm <- obj$data_trans$fit_conc_dose
    }
    if(is.null(log_trans)){
      log_trans <- obj$data_trans$fit_log_conc
    }
    params <- get_params.pkfit(obj)
    #evaluate log-likelihood
    loglike <- log_likelihood(params = params,
                              const_params = NULL,
                              DF = newdata,
                              fit_conc_dose = dose_norm,
                              fit_log_conc = log_trans,
                              force_finite = FALSE,
                              negative = FALSE)
    #calculate BIC
    BIC <- length(params) *
      log(nrow(newdata)) -
      2 * loglike
  }
  return(BIC)
}
