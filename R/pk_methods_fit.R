#' Fit PK model(s) for a `pk` object
#'
#' @param obj A [pk] object.
#' @return The same [pk] object, with element `fit` added to each `$stat_model`
#'   element, reflecting the fitted parameters for the corresponding model.
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


  suppress.messages <- obj$data_settings$suppress.messages
  #For each model:
  for (this_model in names(obj$stat_model)){

    if(obj$stat_model[[this_model]]$status %in% "continue"){

      #Pull par_DF to get which params to optimize, bounds, and starting values
      par_DF <- obj$stat_model[[this_model]]$par_DF
      #Do the same for sigma_DF (the data frame of error SDs)
      sigma_DF <- obj$stat_error_model$sigma_DF

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

      fitdata <- obj$data
      #Now call optimx::optimx() and do the fit
      optimx_out <- tryCatch({

        do.call(
          optimx::optimx,
          args = c(
            #
            list(par = opt_params,
                 fn = log_likelihood,
                 lower = lower_params,
                 upper = upper_params),
            #method and control
            obj$optimx_settings,
            #... additional args to log_likelihood
            list(
              const_params = const_params,
              fitdata = fitdata,
              data_sigma_group = obj$stat_error_model$data_sigma_group,
              modelfun = obj$stat_model[[this_model]]$conc_fun,
              scales_conc = obj$scales$conc,
              negative = TRUE,
              force_finite = TRUE
            ) #end list()
          ) #end args = c()
        ) #end do.call

      },
      error = function(err){
        return(paste0("Error from optimx::optimx(): ",
                      err$message))
      })

      #Save the fitting results for this model
      obj$stat_model[[this_model]]$fit <- optimx_out

    }else{ #if status for this model was "abort", then abort fit and return NULL
      obj$stat_model[[this_model]]$fit <- "Fit aborted because status %in% 'abort'"
    }
  } #end loop over models

  obj$status <- status_fit #fitting complete
  return(obj)
}
