#' Do fitting
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
#' @param n_cores Number of cores used for parallel computing.
#' @param rate_names The names of the rate units. Leave NULL to utilize default 1/hour.
#' @param ... Additional arguments. Not currently in use.
#' @return The same [pk] object, with element `fit` containing the fitted
#'   results for each model in `stat_model`.
#' @export
#' @import multidplyr purrr
#' @author Caroline Ring, Gilberto Padilla Mercado
do_fit.pk <- function(obj, n_cores = NULL, rate_names = NULL, ...){
  #check status
  objname <- deparse(substitute(obj))
  status <- obj$status

  if(status >= status_fit){
    warning(paste0(objname,
                   " current status is ",
                   status,
                   ". do_fit() will reset its status to ",
                   status_fit,
                   ". Any results from later workflow stages will be lost."))
  }

  #if preprocessing not already done, do it
  if (obj$status < status_preprocess) {
    obj <- do_preprocess(obj)
  }

  #if data summary not already done, do it
  if (obj$status < status_data_info) {
    obj <- do_data_info(obj)
  }

  #if prefitting not already done, do it
  if (obj$status < status_prefit) {
    obj <- do_prefit(obj)
  }

  if (!rlang::is_installed("multidplyr")) { n_cores <- NULL }



  #pull the non-excluded observations for fitting
  data <- subset(obj$data, exclude %in% FALSE)
  data_sigma_group <- data$data_sigma_group

  data_group <- get_data_group(obj)
  data_group_vars <- sapply(data_group,
                               rlang::as_label)

  # Make all the other data prior to loading
  par_DF <- obj$prefit$par_DF
  sigma_DF <- obj$prefit$stat_error_model$sigma_DF
  fit_check_DF <- obj$prefit$fit_check


  # nest the necessary data frames...
  data_nest <- data %>%
    tidyr::nest(data = !tidyselect::all_of(data_group_vars))

  par_DF_nest <- par_DF %>%
    tidyr::nest(par_DF = !tidyselect::all_of(c("model", data_group_vars)))

  sigma_DF_nest <-  sigma_DF %>%
    tidyr::nest(sigma_DF = !tidyselect::all_of(data_group_vars))

  fit_check <- fit_check_DF %>%
    dplyr::select(!c(n_par, n_sigma, n_detect, n_par_opt, fit_reason))

  # merge it all together

  info_nest <- suppressMessages(dplyr::inner_join(
    dplyr::inner_join(
      dplyr::inner_join(data_nest,
                        par_DF_nest,
                        by = c(data_group_vars)),
      sigma_DF_nest,
      by = c(data_group_vars)),
    fit_check,
    by = c("model", data_group_vars)) %>%
    dplyr::relocate(model, .after = data_group_vars[-1]) )

  # Set the options for Parallel Computing
  # First condition if it is FALSE don't use parallel computing (takes much longer though)
  #

  if (is.numeric(n_cores)) {
    message(paste0("do_fit.pk(): Trying to divide processes into ", n_cores, " processing cores"))
    total_cores <- parallel::detectCores()
    if (total_cores <= n_cores & total_cores > 1) {
      n_cores  <- total_cores - 1
      message(paste0("do_fit.pk():To ensure other programs & processes are still able to run, ",
                     "n_cores has been set to ", n_cores))
    } else if (total_cores == 1) {
      n_cores = total_cores
    } else {
      n_cores <- n_cores
    }
    message(paste0("do_fit.pk(): ", n_cores, " processing cores allocated."))
    cluster <- multidplyr::new_cluster(n_cores)
    if (any(.packages(all.available = TRUE) %in% "invivoPKfit")) {
      multidplyr::cluster_call(cluster, library(invivoPKfit))
    } else {
      multidplyr::cluster_call(cluster, devtools::load_all())
    }

    multidplyr::cluster_copy(cluster, "obj")

    fit_out <- info_nest %>%
      dplyr::group_by(!!!data_group, model) %>% multidplyr::partition(cluster)
    tidy_fit <- fit_out %>%
      dplyr::summarize(
        fit = purrr::pmap(.l = dplyr::pick(tidyselect::everything()),
                          .f = fit_group,
                          this_model = model,
                          settings_optimx = get_settings_optimx(obj),
                          modelfun = obj$stat_model[[model]]$conc_fun,
                          dose_norm = obj$scales$conc$dose_norm,
                          log10_trans = obj$scales$conc$log10_trans,
                          suppress.messages = TRUE)) %>% dplyr::collect()
  } else {
    fit_out <- info_nest %>%
      dplyr::group_by(!!!data_group, model)

    tidy_fit <- fit_out %>%
      dplyr::summarize(fit = purrr::pmap(.l = dplyr::pick(tidyselect::everything()),
                                         .f = fit_group,
                                         this_model = model,
                                         settings_optimx = get_settings_optimx(obj),
                                         modelfun = obj$stat_model[[model]]$conc_fun,
                                         dose_norm = obj$scales$conc$dose_norm,
                                         log10_trans = obj$scales$conc$log10_trans,
                                         suppress.messages = TRUE))
  }


  # Need to convert rates to perHour
  # Take rate_names
  # Parameter names don't matter, all rates should have consistent param_unit
  message(paste0("do_fit.pk(): Now converting all rate constants to units of 1/hour, ",
                 "in case time has been scaled to units other than hours before fitting"))
  rate_names <- par_DF %>% dplyr::select(!!!obj$data_group,
                                        param_name,
                                        param_units) %>%
    dplyr::filter(stringr::str_detect(param_units, "^1/")) %>%
    dplyr::mutate(Time_trans.Units = stringr::str_remove(param_units, "^1/")) %>%
    dplyr::distinct()
  # Get a simple data_group and conversion rate data frame
  rate_conversion <- rate_names %>%
    dplyr::select(-param_name) %>% # Keep Time_trans.Units for a join
    dplyr::group_by(Time_trans.Units) %>%
    dplyr::mutate(to_perhour = convert_time(1,
                                            from = Time_trans.Units,
                                            to = "hours",
                                            inverse = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()
  # Extract names
  rate_names <- rate_names %>% dplyr::pull(param_name)

  tidy_fit <- tidy_fit %>%
    tidyr::unnest(fit)
  orig_names <- names(tidy_fit) # Saving the original names

  # Convert the rates

  tidy_fit <- suppressMessages(
    tidy_fit %>%
      dplyr::left_join(rate_conversion,
                       by = c(data_group_vars)) %>%
      dplyr::mutate(across(contains(rate_names),
                           \(x) x * to_perhour)))


  # Add optimx bad fits
  obj$conv_not_zero <- tidy_fit %>%
    dplyr::select(-tidyselect::starts_with("sigma")) %>%
    dplyr::filter(convcode != 0)

  # Add parameter fit flags
  obj$params_atBounds <- tidy_fit %>%
    dplyr::select(-(value:xtime), -tidyselect::starts_with("sigma")) %>%
    tidyr::pivot_longer(cols = tidyselect::where(is.numeric),
                        names_to = "param_name",
                        values_to = "param_value") %>%
    dplyr::inner_join(obj$prefit$par_DF) %>% filter(optimize_param) %>%
    mutate(at_bound = ifelse(
      param_value != lower_bound & param_value != upper_bound,
      "Not at bound", ifelse(
        param_value == lower_bound, "LOWER BOUND",
        "UPPER BOUND")
    ))

  obj$fit <- suppressMessages(tidy_fit %>%
                                 dplyr::select(contains(orig_names)) %>%
                                 dplyr::filter(convcode != -9999) %>%
                                 dplyr::group_by(!!!obj$data_group, model) %>%
                                 tidyr::nest(.key = "fit") %>% ungroup())


  obj$status <- status_fit #fitting complete
  message("do_fit.pk: Fitting complete")
  return(obj)
}
