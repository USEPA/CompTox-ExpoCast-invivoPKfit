#' Get parameters to be fitted
#'
#' Given a [pk()] object, get a `data.frame` of parameters and whether they
#' are to be estimated from the data and/or used.
#'
#' @param obj A [pk()] object.
#' @return A `data.frame` with the following variables:
#' - `param_name` Character: the name of each model parameter.
#' - `optimize_param` Logical: Whether the parameter is to be
#'   optimized/estimated from the data (TRUE), or not (FALSE).
#' - `use_param` Logical: Whether the parameter gets used in the model at all
#'   (TRUE) or not (FALSE).
#'
#'   By default, `use_param` and `opt_param` will be the same -- for each
#'   parameter, both `TRUE` or both `FALSE`. However, both variables are
#'   retained to allow for the possibility of holding one or more parameters
#'   constant and optimizing the rest (which would be accomplished by setting
#'   `optimize_param  = FALSE` but keeping `use_param = TRUE` for the parameters
#'   to be held constant).
get_params.pk <- function(obj){
  #get parameter names and units, and determine whether to optimize each of
  #these parameters or not
  par_DF <- do.call(get_opt_params,
                    list("model" = obj$model_info$model,
                         "fitdata" = obj$data,
                         "pool_sigma" = obj$analysis_type %in% "Pooled",
                         "suppress.messages" = suppress.messages))
  return(par_DF)
}

#' Get parameter bounds and starting values
#'
#' @param obj A [pk()] object.

#'@param get_lower_args Any additional arguments to [get_lower_bounds()] (other
#'  than `model` and `fitdata`, which are always passed). Default NULL to accept
#'  the default arguments for [get_lower_bounds()].
#'@param get_upper_args Any additional arguments to [get_upper_bounds()] (other
#'  than `model` and `fitdata`, which are always passed). Default NULL to accept
#' default arguments for [get_upper_bounds()].
#'  @param get_starts_args Any additional arguments to [get_starts()] (other than
#'  `model` and `fitdata`, which are always passed). Default NULL to accept the
#'  default arguments for [get_starts()].
#' @return A `data.frame` the same as `obj$par_DF`, with additional variables `lower_bound`
#'   (numeric, containing the lower bound for each parameter) and
#'   `lower_bound_msg` (character, containing a brief message explaining how the
#'   lower-bound value was calculated).
get_bounds_starts.pk <- function(obj,
                                      get_lower_args = list(Vdist_from_species = FALSE),
                                      get_upper_args = list(Fgutabs_Vdist_from_species = FALSE,
                                                            sigma_from_data = TRUE),
                                      get_starts_args = list(start_from_httk = "all",
                                                             start_from_data = "all")
                                      ){

  for(this_model in names(obj$models)){
  #get lower bounds
  obj[[this_model]]$par_DF <- do.call(get_lower_bounds,
                    c(list("par_DF" = obj[[this_model]]$par_DF,
                           "model" = this_model,
                           "fitdata" = obj$data,
                           "pool_sigma" = obj[[this_model]]$analysis_type %in% "Pooled",
                           "suppress.messages" = suppress.messages),
                      get_lower_args))

  #get upper bounds
  obj[[this_model]]$par_DF <- do.call(get_upper_bounds,
                    c(list("par_DF" = obj[[this_model]]$par_DF,
                           "model" = this_model,
                           "fitdata" = obj$data,
                           "pool_sigma" = obj[[this_model]]$analysis_type %in% "Pooled",
                           "fit_conc_dose" = obj[[this_model]]$data_trans$dose_norm,
                           "fit_log_conc" = obj$data_trans$log10_trans,
                           "suppress.messages" = suppress.messages),
                      get_upper_args)
  )

  #get starting values
  obj$par_DF <- do.call(get_starts,
                    c(list("par_DF" = obj$par_DF,
                           "model" = obj$model_info$model,
                           "fitdata" = obj$data,
                           "pool_sigma" = obj$analysis_type %in% "Pooled",
                           "fit_conc_dose" = obj$data_trans$dose_norm,
                           "fit_log_conc" = obj$data_trans$log10_trans,
                           "suppress.messages" = suppress.messages),
                      get_starts_args))
  }

  return(obj$par_DF)
}

#' Check if there are sufficient detected observations to fit
#'
#' Check whether there are enough detected concentration vs. time observations to fit the
#' requested TK model.
#'
#' @param obj A `pk` object
#' @return A `data.frame` with variables `status`, `n_detect`, `n_opt`. `status`
#'   is the status of the fit setup: either "abort" or "continue". `n_detect` is
#'   the total number of detected observations in the dataset. `n_opt` is the
#'   number of parameters to be optimized. If `n_opt > n_detect`,
#'    then `status =  "abort"`; otherwise, `status = "continue"`.
#' @author Caroline Ring
check_n_detect.pk <- function(obj){
  if("check_n_detect" %in% names(obj)){
    out <- obj$check_n_detect
  }else{
 par_DF <- get_params(obj)
  n_opt <- sum(par_DF$optimize_param %in% TRUE)
  n_detect <- sum(obj$data_info$n_dat$Detect)
  out <- data.frame(n_opt = n_opt,
                           n_detect = n_detect)
  if(n_opt > n_detect){
   out$status <- "abort"
  }else{
  out$status <- "continue"
  }
  }

  return(out)
}

#' Check whether there are sufficient observations in the absorption phase
#'
#' @param obj A `pk` object
#' @return A `data.frame` with variables `fit_kgutabs`, `n_abs_detect`,
#'   `status`.  `fit_kgutabs` is TRUE if the parameter `kgutabs` is to be
#'   optimized; FALSE otherwise. `n_abs_detect` is the total number of detected
#'   observations in the absorption phase (before the empirical time of maximum
#'   dose-normalized concentration).  `status` is the status of the fit setup:
#'   either "abort" or "continue". If `fit_kgutabs %in% TRUE & n_abs_detect < =
#'   2`, then `status =  "abort"`; otherwise, `status = "continue"`.
#' @author Caroline Ring
check_n_abs.pk <- function(obj){
  if("check_n_abs" %in% names(obj)){
    out <- obj$check_n_abs
  }else{
  par_DF <- get_params(obj)
  fit_kgutabs <- "kgutabs" %in% subset(par_DF, optimize_param %in% TRUE)[["param_name"]]
  n_abs <- subset(obj$data_info$n_phase,
                      Phase %in% "absorption")[["Detect"]]
  out <- data.frame(fit_kgutabs = fit_kgutabs,
                    n_abs_detect = n_abs)
 if(fit_kgutabs %in% TRUE &
    n_abs <= 2){
    out$status <- "abort"
  }else{
   out$status <- "continue"
  }
  }

  return(out)
}

#' Check whether there are sufficient observations in the elimination phase
#'
#' @param obj A `pk` object
#' @return A `data.frame` with variables `fit_kelim`, `n_abs_detect`,
#'   `status`.  `fit_kelim` is TRUE if the parameter `kelim` is to be
#'   optimized; FALSE otherwise. `n_elim_detect` is the total number of detected
#'   observations in the absorption phase (before the empirical time of maximum
#'   dose-normalized concentration).  `status` is the status of the fit setup:
#'   either "abort" or "continue". If `fit_kelim %in% TRUE & n_elim_detect < =
#'   2`, then `status =  "abort"`; otherwise, `status = "continue"`.
#' @author Caroline Ring
check_n_elim.pk <- function(obj){
  if("check_n_elim" %in% names(obj)){
    out <- obj$check_n_elim
  }else{
  par_DF <- get_params(obj)
  fit_kelim <- "kelim" %in% subset(par_DF, optimize_param %in% TRUE)[["param_name"]]
  n_elim <- subset(obj$data_info$n_phase,
                  Phase %in% "elimination")[["Detect"]]
  out <- data.frame(fit_kelim = fit_kelim,
                    n_elim_detect = n_elim)
  if(fit_kelim %in% TRUE &
     n_elim <= 2){
    out$status <- "abort"
  }else{
    out$status <- "continue"
  }
  }

  return(out)
}

#' Rescale time
#'
#' If requested, rescale time from hours to some more convenient unit.
#'
#' Time is rescaled based upon the time of the last detected observation in units of hours.
#'
#'  | Range of Last Detect Time in Hours | New Time Units | Range of New Last Detect Time |
#'  | ----------- | -------------- |
#'  | \\[0, 0.5) | minutes | \\[0, 30) |
#'  | \\[0.5, 24) | hours | \\[0.5, 24) |
#'  | \\[24, 720) | days | \\[1, 30) |
#'  | \\[720, 8766) | months | \\[0.99, 12) |
#'  | \\[8766, 87660) | years | \\[1, 10) |
#'  | \\[87660, 876600) | decades | \\[1, 10) |
#'  | \\[876600, Inf) | centuries | \\[1, Inf) |
#'
#'  (At present, no data in CvT requires decades or centuries)
#
#' @param obj A [pk()] object.
#' @return The `data` element of the object, but modified: the original
#'   time (in hours) is now in variable `Time.Hours`; rescaled time is in
#'   `Time`; units of rescaled time are in `Time.Units`.
#' @author Caroline Ring
rescale_time.pk <- function(obj){
  #Rescale time if so requested
  #Save the original time in units of hours
  obj$data$Time.Hours <- obj$data$Time
  #Now do the rescale
  if(obj$data_trans$rescale_time %in% TRUE){
    last_detect_time <- obj$data_info$last_detect_time

    new_time_units <- cut(last_detect_time,
        breaks = c(0,
                   0.5,
                   24,
                   720,
                   8766,
                   87660,
                   876600,
                   Inf),
        labels = c("minutes",
                   "hours",
                   "days",
                   "months",
                   "years",
                   "decades",
                   "centuries"),
        right = FALSE)

    obj$data$Time <- convert_time(x = obj$data$Time.Hours,
                                 from = "hours",
                                 to = new_time_units,
                                 inverse = FALSE)

  }else{
    new_time_units <- "hours"
  }

  obj$data$Time.Units <- new_time_units

  return(obj$data)
}
