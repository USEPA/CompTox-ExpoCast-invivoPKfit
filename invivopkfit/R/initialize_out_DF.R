initialize_out_DF <- function(par_DF,
                              msg = "",
                              fitted_types = c("Fitted mean",
                                               "Fitted std dev"),
                              studies_analyzed,
                              refs_analyzed,
                              analysis_type,
                              fit_conc_dose,
                              fit_log_conc,
                              n_routes,
                              n_media,
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

  out_DF[, fitted_types] <- NA_real_

  out_DF$Studies.Analyzed <- studies_analyzed
  out_DF$References.Analyzed <- refs_analyzed
  out_DF$Data.Analyzed <- analysis_type

  out_DF$fit_conc_dose <- fit_conc_dose
  out_DF$fit_log_conc <- fit_log_conc

  #fill in the loglike and AIC with NA s since no fit was done
  out_DF$LogLikelihood <-  NA_real_
  out_DF$AIC <-  NA_real_


  out_DF$flag <- NA_character_
  #Record the unique routes in this dataset
  #Route info provides context for why some parameters were/were not estimated
  out_DF$N_Routes <- n_routes
  #Record the unique media in this dataset
  out_DF$N_Media <- n_media

  out_DF$message <- msg
  out_DF$fevals <- NA_integer_
  out_DF$convcode <- NA_integer_
  out_DF$niter <- NA_integer_
  out_DF$method <- optimx_args$method
  #record control params
  out_DF[paste0("control_",
                names(optimx_args$control))] <- optimx_args$control

  if(!suppress.messages){
    message(msg)
  }
  return(out_DF)
}
