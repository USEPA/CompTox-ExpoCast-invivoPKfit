
#' Create a `pksetup` object
#'
#' @param pkdata An object of class `pkdata`
#' @param analysis_type Which type of analysis to do. One of "Joint",
#'   "Separate", or "Pooled" (default "Joint"). "Joint" will fit data from
#'   multiple studies to one single concentration vs. time model, allowing a
#'   different residual error variance for each study. "Separate" will fit data
#'   from multiple studies to multiple concentration vs. time models, one for
#'   each study. "Pooled" will fit data from multiple studies to one single
#'   concentration vs. time model, assuming the same residual error variance
#'   applies to all studies.
#' @param model_info Which toxicokinetic (concentration vs. time) model to fit. A
#'   named list with two elements. The first element is "model", which may be
#'   one of "flat", "1compartment", or "2compartment". Details on each model can
#'   be found in the documentation for [cp_flat()], [cp_1compartment()], and
#'   [cp_2compartment()], respectively. The second element is "model.type",
#'   which may be one of "analytic" (default) or "full". "analytic" will use the
#'   analytic solution for the flat, 1-compartment, or 2-compartment models.
#'   "full" will use the full ODE model (solved numerically) and is currently
#'   only implemented for `model = "1compartment"`. `model.type = 'analytic'` is
#'   strongly recommended; the solutions are exact and will be much faster.
#' @param data_trans Whether and how the data should be transformed before
#'   fitting. A named list with elements "dose_norm", "log10_trans",
#'   "rescale_time". Default is `list(dose_norm = TRUE, log10_trans = FALSE,
#'   rescale_time = TRUE)`. "dose_norm" controls whether concentrations should
#'   be normalized to dose before fitting. "log10_trans" controls whether
#'   concentrations should be log10-transformed before fitting. "rescale_time"
#'   controls whether time (in hours) should be automatically rescaled before
#'   fitting (to units determined by the maximum value of time).
#' @param optim_control Arguments for the optimizer ([optimx::optimx()]). A
#'   named list whose elements correspond to arguments for [optimx::optimx()]
#'   (other than "par", "fn", "gr", "hess", "lower", or "upper", which will all
#'   be assigned automatically). Default is `list(method = "bobyqa", itnmax =
#'   1e6, control = list(kkt = FALSE))`.
#' @return An object of class `pksetup`.
#' @export
#' @author Caroline Ring
pksetup <- function(pkdata = NULL,
                    analysis_type = "Joint",
                    model_info = list(model = "flat",
                                      model.type = "analytic"),
                    data_trans = list(dose_norm = TRUE,
                                      log10_trans = FALSE,
                                      rescale_time = TRUE),
                    optim_control = list(
                      method = "bobyqa",
                      itnmax = 1e6,
                      control = list(kkt = FALSE)
                    ),
                    get_lower_args = list(Vdist_from_species = FALSE),
                    get_upper_args = list(Fgutabs_Vdist_from_species = FALSE,
                                          sigma_from_data = TRUE),
                    get_starts_args = list(start_from_httk = "all",
                                           start_from_data = "all")
){

  if(is.null(pkdata)){
    #create an empty data set
    pkdata <- pkdata()
  }

  obj <- c(pkdata,
                   list("analysis_type" = analysis_type,
                        "model_info" = model_info,
                        "data_trans" = data_trans,
                        "optim_control" = optim_control))

  class(obj) <- c(class(obj),
                                "pksetup")

  #now add:
  #parameters
  obj$par_DF <- get_params.pksetup(obj)
  #checks
  obj$check_n_detect <- check_n_detect.pksetup(obj)
  obj$check_n_abs <- check_n_abs.pksetup(obj)
  obj$check_n_elim <- check_n_elim.pksetup(obj)
  #overall fit status and reason
  obj$status <- ifelse(
    obj$check_n_detect$status %in% "abort",
    "abort",
    ifelse(
      obj$check_n_abs$status %in% "abort",
      "abort",
      ifelse(
        obj$check_n_elim$status %in% "abort",
        "abort",
        "continue"
      )
    )
  )

  obj$status_reason <- ifelse(
    obj$check_n_detect$status %in% "abort",
    "Fewer detected observations than parameters to optimize",
    ifelse(
      obj$check_n_abs$status %in% "abort",
      "Fewer than 3 detected observations during absorption phase; cannot optimize kgutabs",
      ifelse(
        obj$check_n_elim$status %in% "abort",
        "Fewer than 3 detected observations during elimination phase; cannot optimize kelim",
        "Sufficient observations to optimize all requested parameters"
      )
    )
  )

  #rescale time
  obj$data <- rescale_time.pksetup(obj)

  #get bounds & starting values for parameter optimization
 obj$par_DF <- get_bounds_starts.pksetup(obj)


  return(obj)
}

#' Get parameters to be fitted
#'
#' Given a [pksetup()] object, get a `data.frame` of parameters and whether they
#' are to be estimated from the data and/or used.
#'
#' @param obj A [pksetup()] object.
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
get_params.pksetup <- function(obj){
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
#' @param obj A [pksetup()] object.

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
get_bounds_starts.pksetup <- function(obj,
                                      get_lower_args = list(Vdist_from_species = FALSE),
                                      get_upper_args = list(Fgutabs_Vdist_from_species = FALSE,
                                                            sigma_from_data = TRUE),
                                      get_starts_args = list(start_from_httk = "all",
                                                             start_from_data = "all")
                                      ){
  #get lower bounds
  obj$par_DF <- do.call(get_lower_bounds,
                    c(list("par_DF" = obj$par_DF,
                           "model" = obj$model_info$model,
                           "fitdata" = obj$data,
                           "pool_sigma" = obj$analysis_type %in% "Pooled",
                           "suppress.messages" = suppress.messages),
                      get_lower_args))

  #get upper bounds
  obj$par_DF <- do.call(get_upper_bounds,
                    c(list("par_DF" = obj$par_DF,
                           "model" = obj$model_info$model,
                           "fitdata" = obj$data,
                           "pool_sigma" = obj$analysis_type %in% "Pooled",
                           "fit_conc_dose" = obj$data_trans$dose_norm,
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

  return(obj$par_DF)
}

#' Check if there are sufficient detected observations to fit
#'
#' Check whether there are enough detected concentration vs. time observations to fit the
#' requested TK model.
#'
#' @param obj A `pksetup` object
#' @return A `data.frame` with variables `status`, `n_detect`, `n_opt`. `status`
#'   is the status of the fit setup: either "abort" or "continue". `n_detect` is
#'   the total number of detected observations in the dataset. `n_opt` is the
#'   number of parameters to be optimized. If `n_opt > n_detect`,
#'    then `status =  "abort"`; otherwise, `status = "continue"`.
#' @author Caroline Ring
check_n_detect.pksetup <- function(obj){
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

  return(out)
}

#' Check whether there are sufficient observations in the absorption phase
#'
#' @param obj A `pksetup` object
#' @return A `data.frame` with variables `fit_kgutabs`, `n_abs_detect`,
#'   `status`.  `fit_kgutabs` is TRUE if the parameter `kgutabs` is to be
#'   optimized; FALSE otherwise. `n_abs_detect` is the total number of detected
#'   observations in the absorption phase (before the empirical time of maximum
#'   dose-normalized concentration).  `status` is the status of the fit setup:
#'   either "abort" or "continue". If `fit_kgutabs %in% TRUE & n_abs_detect < =
#'   2`, then `status =  "abort"`; otherwise, `status = "continue"`.
#' @author Caroline Ring
check_n_abs.pksetup <- function(obj){
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

  return(out)
}

#' Check whether there are sufficient observations in the elimination phase
#'
#' @param obj A `pksetup` object
#' @return A `data.frame` with variables `fit_kelim`, `n_abs_detect`,
#'   `status`.  `fit_kelim` is TRUE if the parameter `kelim` is to be
#'   optimized; FALSE otherwise. `n_elim_detect` is the total number of detected
#'   observations in the absorption phase (before the empirical time of maximum
#'   dose-normalized concentration).  `status` is the status of the fit setup:
#'   either "abort" or "continue". If `fit_kelim %in% TRUE & n_elim_detect < =
#'   2`, then `status =  "abort"`; otherwise, `status = "continue"`.
#' @author Caroline Ring
check_n_elim.pksetup <- function(obj){
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
#' @param obj A [pksetup()] object.
#' @return The same object, but with the `data` element modified: the original
#'   time (in hours) is now in variable `Time.Hours`; rescaled time is in
#'   `Time`; units of rescaled time are in `Time.Units`.
#' @author Caroline Ring
rescale_time.pksetup <- function(obj){
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
