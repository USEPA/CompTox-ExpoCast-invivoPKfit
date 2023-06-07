#' Fitting
#'
#' Fit PK model(s) for a `pk` object
#'
#' This function estimates the parameters for each model in `stat_model` from
#' the data, using numerical optimization implemented in [optimx::optimx()]. The
#' optimization is done by maximizing the log-likelihood function implemented in
#' [log_likelihood()] (technically, by minimizing the negative log-likelihood).
#' Only the non-excluded observations are used.
#'
#' Due to limitations of [optimx::optimx()], the log-likelihood function is
#' forced to return finite values during this optimization. Impossible
#' combinations of parameters (e.g., parameter values that produce negative
#' predicted concentrations) should have a log-likelihood of `-Inf`, but due to
#' this limitation, they instead have a log-likelihood of `-Machine.doublexmax`.
#' This limitation means that the log-likelihood function is flat in regions of
#' impossible parameter values. It is unlikely, but possible, that the optimizer
#' might get "stuck" in such a flat region -- report convergence, but return a
#' "bad" set of parameter values that produces non-physical predictions.
#'
#' Before trusting the results of any fit, it is recommended to check the
#' log-likelihood using [logLik()] and the Akaike Information Criterion using
#' [AIC()], which check the log-likelihood *without* forcing it to return finite
#' values.
#'
#' @param obj A [pk] object.
#' @return The same [pk] object, with element `fit` containing the fitted
#'   results for each model in `stat_model`.
#' @export
#' @author Caroline Ring
fit.pk <- function(obj){
  #check status
  objname <- deparse(substitute(obj))
  status <- obj$status
  if(status >= status_fit){
    warning(paste0(objname,
                   " current status is ",
                   status,
                   ". fit() will reset its status to ",
                   status_fit,
                   ". Any results from later workflow stages will be lost."))
  }

  #if preprocessing not already done, do it
  if(obj$status < status_preprocess){
    obj <- preprocess_data(obj)
  }

  if(obj$status < status_data_info){
    obj <- data_info(obj)
  }

  if(obj$status < status_prefit){
    obj <- prefit(obj)
  }

  #pull the non-excluded observations for fitting
  data <- subset(obj$data, exclude %in% FALSE)
  #pull the sigma parameter corresponding to each non-excluded observation
  data_sigma_group <- obj$prefit$stat_error_model$data_sigma_group[obj$data$exclude %in% FALSE]

  suppress.messages <- obj$settings_preprocess$suppress.messages
  #For each model:
  for (this_model in names(obj$stat_model)){

    if(obj$prefit[[this_model]]$fit_decision %in% "continue"){

      #Pull par_DF to get which params to optimize, bounds, and starting values
      par_DF <- obj$prefit[[this_model]]$par_DF
      #Do the same for sigma_DF (the data frame of error SDs)
      sigma_DF <- obj$prefit$stat_error_model$sigma_DF

      #Rowbind par_DF and sigma_DF
      par_DF <- rbind(par_DF,
                      sigma_DF)

      #get params to be optimized with their starting points
      opt_params <- par_DF[par_DF$optimize_param %in% TRUE,
                           "start"]
      names(opt_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                  "param_name"]

      #param lower bounds (only on params to be optimized)
      lower_params <- par_DF[par_DF$optimize_param %in% TRUE,
                             "lower_bound"]
      names(lower_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                    "param_name"]

      #param upper bounds (only on params to be optimized)
      upper_params <- par_DF[par_DF$optimize_param %in% TRUE,
                             "upper_bound"]
      names(upper_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                    "param_name"]

      #get params to be held constant, if any
      if(any(par_DF$optimize_param %in% FALSE &
             par_DF$use_param %in% TRUE)){
        const_params <- par_DF[par_DF$optimize_param %in% FALSE &
                                 par_DF$use_param %in% TRUE,
                               "start"]

        names(const_params) <- par_DF[par_DF$optimize_param %in% FALSE &
                                        par_DF$use_param %in% TRUE,
                                      "param_name"]
      }else{
        const_params <- NULL
      }


      if(suppress.messages %in% FALSE){
        message(paste("fit.pk(): Fitting model",
                this_model,
                "using optimx::optimx()"))
      }

      #Now call optimx::optimx() and do the fit
      suppressWarnings(
        optimx_out <- tryCatch({

        tmp <- do.call(
          optimx::optimx,
          args = c(
            #
            list(par = opt_params,
                 fn = log_likelihood,
                 lower = lower_params,
                 upper = upper_params),
            #method and control
            obj$settings_optimx,
            #... additional args to log_likelihood
            list(
              const_params = const_params,
              data = data,
              data_sigma_group = data_sigma_group,
              modelfun = obj$stat_model[[this_model]]$conc_fun,
              dose_norm = obj$scales$conc$dose_norm,
              log10_trans = obj$scales$conc$log10_trans,
              negative = TRUE,
              force_finite = TRUE,
              suppress.messages = suppress.messages
            ) #end list()
          ) #end args = c()
        ) #end do.call

        # tmp <- cbind(tmp,
        #              attr(tmp, "details"))
        # #in case convcode is 9999, this implies an error in calling stat::optim
        # tmp[tmp$convcode %in% 9999, "message"] <- "Error from stat::optim() called by optimx::optimx()"

        tmp

      },
      error = function(err){
        return(paste0("Error from optimx::optimx(): ",
                      err$message))
      }) #end tryCatch()
      ) #end suppressWarnings()

      #Save the fitting results for this model
      obj$fit[[this_model]] <- optimx_out

    }else{ #if status for this model was "abort", then abort fit and return NULL
      obj$fit[[this_model]] <- "Fit aborted because status %in% 'abort'"
    }
  } #end loop over models

  obj$status <- status_fit #fitting complete
  return(obj)
}
