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
#' @return The same [pk] object, with element `fit` containing the fitted
#'   results for each model in `stat_model`.
#' @export
#' @author Caroline Ring
do_fit.pk <- function(obj, n_cores = NULL, rate_names = NULL){
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


  #nest the necessary data frames...
  data_nest <- data %>%
    tidyr::nest(data = !tidyselect::all_of(data_group_vars))

  par_DF_nest <- par_DF %>%
    tidyr::nest(par_DF = !tidyselect::all_of(c("model", data_group_vars)))

  sigma_DF_nest <-  sigma_DF %>%
    tidyr::nest(sigma_DF = !tidyselect::all_of(data_group_vars))

  fit_check <- fit_check_DF %>%
    dplyr::select(!c(n_par, n_sigma, n_detect, n_par_opt, fit_reason))

  #merge it all together

  info_nest <- dplyr::inner_join(
    dplyr::inner_join(
      dplyr::inner_join(data_nest,
                        par_DF_nest,
                        by = c(data_group_vars)),
      sigma_DF_nest,
      by = c(data_group_vars)),
    fit_check,
    by = c("model", data_group_vars)) %>%
    dplyr::relocate(model, .after = data_group_vars[-1])

  total_cores <- parallel::detectCores()
  if (is.null(n_cores)) {
    if (total_cores > 4) {
      n_cores <- total_cores - 2
    } else if (total_cores >= 2 & total_cores < 4) {
      n_cores <- total_cores - 1
    } else {
      n_cores <- 1
    }
  } else {
    if (total_cores <= n_cores & total_cores > 1) {
      n_cores = total_cores - 1
      message(paste0("To ensure other programs & processes are still to run, ",
              "n_cores has been set to ", n_cores))
    } else if (total_cores == 1) {
      n_cores = total_cores
    } else {
      n_cores = n_cores
    }
  }

  cluster <- multidplyr::new_cluster(n_cores)
  if (any(.packages(all.available = TRUE) %in% "invivoPKfit")) {
    multidplyr::cluster_call(cluster, library(invivoPKfit))
  } else {
    multidplyr::cluster_call(cluster, devtools::load_all())
  }

  multidplyr::cluster_copy(cluster, "obj")

  fit_out <- info_nest %>%
    dplyr::group_by(!!!data_group, model) %>% multidplyr::partition(cluster)
  tidy_fit <- fit_out %>% dplyr::summarize(fit = purrr::pmap(.l = dplyr::pick(tidyselect::everything()),
                                                              .f = fit_group,
                                                              this_model = model,
                                                              settings_optimx = get_settings_optimx(obj),
                                                              modelfun = obj$stat_model[[model]]$conc_fun,
                                                              dose_norm = obj$scales$conc$dose_norm,
                                                              log10_trans = obj$scales$conc$log10_trans,
                                                              suppress.messages = TRUE)) %>%
    dplyr::collect()


  # Need to convert rates to perHour
  # Take rate_names
  # Parameter names don't matter, all rates should have consistent param_unit
  message("Now doing any rate conversions!")
  rate_names <- par_DF %>% dplyr::select(!!!obj$data_group,
                                        param_name,
                                        param_units) %>%
    dplyr::filter(str_detect(param_units, "^1/")) %>%
    dplyr::mutate(Time_trans.Units = str_remove(param_units, "^1/")) %>%
    dplyr::distinct()
  # Get a simple data_group and conversion rate data frame
  rate_conversion <- rate_names %>%
    dplyr::select(-param_name) %>% # Keep Time_trans.Units for a join
    dplyr::mutate(to_perhour = convert_time(1,
                                            from = Time_trans.Units,
                                            to = "hours",
                                            inverse = TRUE)) %>%
    dplyr::distinct()
  # Extract names
  rate_names <- rate_names %>% dplyr::pull(param_name)

  tidy_fit <- tidy_fit %>%
    tidyr::unnest(fit)
  orig_names <- names(tidy_fit) # Saving the original names

  # Convert the rates
  tidy_fit <- tidy_fit %>%
    dplyr::left_join(rate_conversion,
              by = c(data_group_vars)) %>%
    dplyr::mutate(across(contains(rate_names),
                         \(x) x * to_perhour)) %>%
    dplyr::select(contains(orig_names)) %>%
    dplyr::group_by(!!!obj$data_group, model) %>%
    tidyr::nest(.key = "fit")


  obj$fit <- tidy_fit %>%
    dplyr::ungroup()

  obj$status <- status_fit #fitting complete
  return(obj)
}
