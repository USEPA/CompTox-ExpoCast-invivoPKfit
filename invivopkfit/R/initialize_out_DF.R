initialize_out_DF <- function(par_DF,
                              all_data_fit = NULL,
                              msg = "",
                              fitted_types = c("Fitted mean",
                                               "Fitted std dev"),
                              studies_analyzed,
                              refs_analyzed,
                              analysis_type,
                              fit_conc_dose,
                              fit_log_conc,
                              rescale_time,
                              n_routes,
                              n_media,
                              tmax = NA_real_,
                              n_abs = NA_integer_,
                              n_elim = NA_integer_,
                              optimx_args,
                              suppress.messages = TRUE){
  out_DF <- par_DF[par_DF$optimize_param %in% TRUE, ]

  #add NA for lower, upper bounds and start values
  out_DF[, c("lower_bound",
             "lower_bound_msg",
             "upper_bound",
             "upper_bound_msg",
             "start_value",
             "start_value_msg")] <- list(NA_real_,
                                         NA_character_,
                                         NA_real_,
                                         NA_character_,
                                         NA_real_,
                                         NA_character_)

  out_DF$time_units_fitted <- "hours"

  if(is.null(all_data_fit)){
  out_DF[, fitted_types] <- NA_real_
  }else{
    #if we have values for parameters from all_data_fit, keep them
    out_DF[["Fitted mean"]] <- unlist(all_data_fit[out_DF$param_name])
    #but SDs will be NA
    out_DF[["Fitted std dev"]] <- NA_real_
  }

  out_DF$Studies.Analyzed <- studies_analyzed
  out_DF$References.Analyzed <- refs_analyzed
  out_DF$Data.Analyzed <- analysis_type

  out_DF$fit_conc_dose <- fit_conc_dose
  out_DF$fit_log_conc <- fit_log_conc
  out_DF$rescale_time <- rescale_time

  #fill in the loglike and AIC with NA s since no fit was done
  out_DF$LogLikelihood <-  NA_real_
  out_DF$AIC <-  NA_real_


  out_DF$flag <- NA_character_
  #Record the unique routes in this dataset
  #Route info provides context for why some parameters were/were not estimated
  out_DF$N_Routes <- n_routes
  #Record the unique media in this dataset
  out_DF$N_Media <- n_media

  #record tmax, n_abs, n_elim
  out_DF$tmax_oral <- tmax
  out_DF$n_abs <- n_abs
  out_DF$n_elim <- n_elim

  out_DF$message <- msg
  if(is.null(all_data_fit)){
  out_DF$fevals <- NA_real_
  out_DF$convcode <- NA_real_
  out_DF$niter <- NA_real_
  }else{
    out_DF$fevals <- all_data_fit$fevals
    out_DF$convcode <- all_data_fit$convcode
    out_DF$niter <- all_data_fit$niter
  }
  out_DF$method <- optimx_args$method
  #record control params
  out_DF[paste0("control_",
                names(optimx_args$control))] <- optimx_args$control

  if(!suppress.messages){
    message(msg)
  }
  return(out_DF)
}
