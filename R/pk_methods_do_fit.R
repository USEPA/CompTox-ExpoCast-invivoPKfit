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
#' @inheritParams do_preprocess.pk
#' @param async Logical (Default: FALSE). Whether to use for parallel computing with
#'  the [mirai::mirai_map()] framework. If FALSE and [mirai::daemons()] are set, will use
#'  new parallel processing in [purrr::in_parallel()] or if none are set, then sequential
#'  iteration will occur. Progress bars are displayed if TRUE, or if not using parallel
#'  processing.
#' @param rate_names The names of the rate units. Leave NULL to utilize default 1/hour.
#' @return The same [pk] object, with element `fit` containing the fitted
#'   results for each model in `stat_model`.
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
do_fit.pk <- function(obj,
                      async = FALSE,
                      rate_names = NULL,
                      ...) {
  # check status
  objname <- deparse(substitute(obj))
  status <- obj$status

  if (status >= status_fit) {
    warning(
      objname,
      " current status is ",
      status,
      ". do_fit() will reset its status to ",
      status_fit,
      ". Any results from later workflow stages will be lost."
    )
  }
  cli::cli_par()

  # if preprocessing not already done, do it
  if (obj$status < status_preprocess) {
    obj <- do_preprocess(obj)
  }

  # if data summary not already done, do it
  if (obj$status < status_data_info) {
    obj <- do_data_info(obj)
  }

  # if prefitting not already done, do it
  if (obj$status < status_prefit) {
    obj <- do_prefit(obj)
  }

  # Check async capability
  if (async && isFALSE(mirai::daemons_set())) {
    cli::cli_warn(c(
      "!" = paste(
        "Remember to run",
        "{.run [mirai::everywhere(library(invivoPKfit))](mirai::everywhere())} ",
        "after setting number of daemons."
      )
    ))
    mirai::require_daemons() # Will error out of function.
  }

  if (isTRUE(mirai::daemons_set())) {
    cli::cli_inform("Using parallel or asynchronous processing!")
  }

  # pull the non-excluded observations for fitting
  data <- subset(get_data.pk(obj), exclude == FALSE)
  data_sigma_group <- data$data_sigma_group

  data_grp <- get_data_group.pk(obj)
  data_grp_vars <- get_data_group.pk(obj, as_character = TRUE)

  error_group_vars <- get_error_group.pk(obj, as_character = TRUE)

  req_data_vars <- c(
    "DATA_GROUP_ID", # This will help nest all data
    "Route", "Media", "Dose",
    "Conc", "Conc_trans", "Conc_SD",
    "LOQ", "exclude", "Detect", "pLOQ",
    "N_Subjects", "data_sigma_group", "Time_trans"
  )

  # Make all the other data prior to loading
  par_DF <- obj$prefit$par_DF
  sigma_DF <- obj$prefit$sigma_DF
  fit_check_DF <- obj$prefit$fit_check |>
    dplyr::select(!c(n_par, n_sigma, n_detect, n_par_opt, fit_reason))


  # Subset data so I only have the columns needed for fitting
  data <- data[union(req_data_vars, error_group_vars)]
  data_group_keys <- data |>
    dplyr::distinct(DATA_GROUP_ID, !!!data_grp)

  # nest the necessary data frames...
  data_nest <- data |>
    dplyr::group_by(DATA_GROUP_ID) |>
    tidyr::nest(.key = "data") |>
    dplyr::arrange(DATA_GROUP_ID)

  par_DF_nest <- data_group_keys |>
    dplyr::left_join(par_DF, by = c("DATA_GROUP_ID", data_grp_vars)) |>
    dplyr::group_by(model, DATA_GROUP_ID) |>
    tidyr::nest(.key = "par_DF") |>
    dplyr::arrange(DATA_GROUP_ID)

  sigma_DF_nest <- data_group_keys |>
    dplyr::left_join(sigma_DF, by = data_grp_vars) |>
    dplyr::group_by(DATA_GROUP_ID) |>
    tidyr::nest(.key = "sigma_DF") |>
    dplyr::arrange(DATA_GROUP_ID)

  fit_check <- data_group_keys |>
    dplyr::left_join(fit_check_DF, by = data_grp_vars) |>
    dplyr::select(-dplyr::all_of(data_grp_vars)) |>
    dplyr::arrange(DATA_GROUP_ID)

  # make a join-able data.frame with all the possible models
  fun_models <- get_stat_model(obj)

  # merge it all together
  # fit_check > data_nest > par_DF_nest > sigma_DF_nest
  info_nest <- fit_check |>
    dplyr::inner_join(
      data_nest,
      by = "DATA_GROUP_ID",
      relationship = "many-to-one"
    ) |>
    dplyr::inner_join(
      par_DF_nest,
      by = c("DATA_GROUP_ID", "model"),
      relationship = "one-to-one"
    ) |>
    dplyr::inner_join(
      sigma_DF_nest,
      by = "DATA_GROUP_ID",
      relationship = "many-to-one"
    ) |>
    dplyr::inner_join(
      fun_models,
      by = join_by(model),
      relationship = "many-to-one"
    )

  this_settings_optimx <- get_settings_optimx(obj)
  these_methods <- rlang::eval_tidy(this_settings_optimx$method)
  dose_norm <- obj$scales$conc$dose_norm
  log10_trans <- obj$scales$conc$log10_trans

  cli_inform("do_fit.pk(): Begin fitting for model{?s}: {fun_models$model}")

  if (async && isTRUE(mirai::daemons_set())) {
    info_prep <- info_nest |>
      dplyr::mutate(
        settings_optimx = list(this_settings_optimx),
        dose_norm = dose_norm,
        log10_trans = log10_trans,
        suppress.messages = TRUE
      ) |>
      dplyr::select(
        data,
        par_DF,
        sigma_DF,
        fit_decision,
        this_model = model,
        settings_optimx,
        modelfun,
        dose_norm,
        log10_trans,
        suppress.messages
      )

    tidy_fit <- mirai::mirai_map(info_prep, fit_group)[.progress, .stop]
  } else {
    # Things are being passed down as list???
    tidy_fit <- purrr::pmap(
      list(
        data = info_nest$data,
        par_DF = info_nest$par_DF,
        sigma_DF = info_nest$sigma_DF,
        fit_decision = info_nest$fit_decision,
        this_model = info_nest$model,
        settings_optimx = list(this_settings_optimx),
        modelfun = info_nest$modelfun,
        dose_norm = dose_norm,
        log10_trans = log10_trans,
        suppress.messages = TRUE
      ),
      .f = purrr::in_parallel(
        \(data, par_DF, sigma_DF,
          fit_decision, this_model,
          settings_optimx, modelfun,
          dose_norm, log10_trans,
          suppress.messages) {
          fit_group(
            data, par_DF, sigma_DF,
            fit_decision, this_model,
            settings_optimx, modelfun,
            dose_norm, log10_trans,
            suppress.messages
          )
        },
        fit_group = fit_group
      ),
      .progress = list(
        type = "iterator",
        format = "Fitting... {cli::pb_bar} {cli::pb_current}/{cli::pb_total} [{cli::pb_elapsed}]",
        clear = FALSE
      )
    )
  }
  # Need to cbind this list with the original data
  tidy_fit <- info_nest |>
    dplyr::mutate(fit = tidy_fit)

  tidy_fit <- tidy_fit |>
    dplyr::mutate(
      fit = purrr::map(
        fit,
        function(x) {
          x |>
            dplyr::mutate(message = as.character(message)) |>
            tidyr::pivot_longer(
              cols = !c(method, value:message),
              names_to = "param_name",
              values_to = "estimate"
            ) |>
            dplyr::relocate(
              param_name, estimate, convergence, value,
              .after = method
            ) |>
            dplyr::mutate(estimate = as.list(estimate))
        }
      )
    )

  # Unnest
  tidy_fit <- tidy_fit |>
    tidyr::unnest(fit) |>
    dplyr::select(DATA_GROUP_ID, model, method:message)

  tidy_sigmas <- sigma_DF_nest |>
    tidyr::unnest(sigma_DF) |>
    dplyr::select(DATA_GROUP_ID, !!!data_grp, param_name:start) |>
    dplyr::mutate(start = as.list(start)) |>
    dplyr::cross_join(
      expand.grid(list(model = unique(par_DF$model), method = these_methods)),
      copy = TRUE
    ) |>
    dplyr::relocate(model, method, .after = DATA_GROUP_ID)

  tidy_params <- par_DF_nest |>
    tidyr::unnest(par_DF) |>
    dplyr::cross_join(
      list(method = these_methods),
      copy = TRUE
    ) |>
    dplyr::relocate(model, method, .after = DATA_GROUP_ID)

  tidy_fit <- dplyr::bind_rows(tidy_params, tidy_sigmas) |>
    dplyr::full_join(tidy_fit, by = c("DATA_GROUP_ID", "model", "method", "param_name")) |>
    dplyr::filter(use_param) |>
    dplyr::group_by(DATA_GROUP_ID, model, method) |>
    tidyr::fill(convergence:message, .direction = "downup") |>
    dplyr::ungroup() |>
    # because start is a list, need toedo type conversion
    dplyr::mutate(
      estimate = ifelse(
        optimize_param == FALSE & use_param == TRUE & !convergence %in% c(9999, -9999),
        start, estimate
      )
    )

  # Take rate_names
  # Parameter names don't matter, all rates should have consistent param_unit
  cli_inform(paste(
    "do_fit.pk(): Now converting all rate constant estimates to units of 1/hour,",
    "in case time has been scaled to units other than hours before fitting."
  ))

  rate_names <- par_DF |>
    dplyr::distinct(
      !!!data_grp,
      model,
      param_name,
      param_units
    ) |>
    dplyr::filter(startsWith(param_units, "1/")) |>
    dplyr::mutate(Time_trans.Units = stringr::str_remove(param_units, "^1/"))

  # Get a simple data_group and conversion rate data frame
  rate_conversion <- rate_names |>
    dplyr::group_by(Time_trans.Units) |>
    dplyr::mutate(to_perhour = convert_time(1,
      from = Time_trans.Units,
      to = "hours",
      inverse = TRUE
    )) |>
    dplyr::ungroup() |>
    dplyr::distinct()

  # Convert the rates
  tidy_fit <- tidy_fit |>
    dplyr::left_join(rate_conversion, by = c(
      data_grp_vars,
      "model",
      "param_name",
      "param_units"
    ))

  rate_param_rows <- with(tidy_fit, which(!is.na(to_perhour)))

  tidy_fit[rate_param_rows, ][["start"]] <- as.list(
    unlist(tidy_fit[rate_param_rows, ][["start"]]) * tidy_fit[rate_param_rows, ][["to_perhour"]]
  )
  tidy_fit[rate_param_rows, ][["estimate"]] <- as.list(
    unlist(tidy_fit[rate_param_rows, ][["estimate"]]) * tidy_fit[rate_param_rows, ][["to_perhour"]]
  )
  tidy_fit <- dplyr::select(tidy_fit, -to_perhour)

  # Changing the final fit form so that everything is similar par_DF

  # Add parameter fit flags
  obj$fit <- tidy_fit |>
    dplyr::mutate(
      at_bound = dplyr::case_when(
        identical(estimate, lower_bound) ~ "AT LOWER BOUND",
        identical(estimate, upper_bound) ~ "AT UPPER BOUND",
        identical(estimate, start) ~ "AT START",
        .default = "Not at bound"
      )
    ) |>
    dplyr::select(-c(lower_bound, upper_bound)) |>
    dplyr::distinct()

  obj$status <- status_fit # fitting complete
  cli_inform(c("v" = "do_fit.pk(): Fitting complete!"))
  return(obj)
}
