#' Fit a single group of data
#'
#' @param data A single group of data
#' @param sigma_DF sigma_DF for a single group of data
#' @param par_DF par_DF for a single group of data
#' @param this_model Name of the `pk_model` object to fit
#' @param fit_decision Whether the fit is able to be calculated or excluded.
#' @param settings_optimx The settings for optimization.
#' @param modelfun Name of the model concentration function
#' @param dose_norm TRUE or FALSE -- whether to dose-normalize concentrations
#'   before evaluating log-likelihood
#' @param log10_trans TRUE or FALSE -- whether to 1og10-transform concentrations
#'   before evaluating log-likelihood
#' @param suppress.messages TRUE or FALSE -- whether to suppress messages or emit them
#' @return An object of class `optimx` (i.e. a data.frame with fit results)
fit_group <- function(data,
                      par_DF,
                      sigma_DF,
                      fit_decision,
                      this_model,
                      settings_optimx,
                      modelfun,
                      dose_norm,
                      log10_trans,
                      suppress.messages) {
  # get params to be held constant, if any
  if (any(par_DF$optimize_param %in% FALSE & par_DF$use_param %in% TRUE)) {
    const_params <- par_DF[which(!par_DF$optimize_param & par_DF$use_param),
                           "start", drop = TRUE]
    names(const_params) <- par_DF[
      which(!par_DF$optimize_param & par_DF$use_param),
      "param_name",
      drop = TRUE
    ]
  } else {
    const_params <- NULL
  }

  # Now par_DF should only be optimized parameters, and unnest the start values
  par_DF <- par_DF |> dplyr::filter(optimize_param) |>
    tidyr::unnest(cols = start)

  # Rowbind par_DF and sigma_DF
  par_DF <- dplyr::bind_rows(par_DF, sigma_DF)

  # get params to be optimized with their starting points
  # since this series of assignments only return ordered vectors, possible to optimize a little
  # also optimized parameter names was repeated, simply put into variable
  opt_param_names <- par_DF[which(par_DF$optimize_param), "param_name",
                            drop = TRUE]

  opt_params <- par_DF[which(par_DF$optimize_param), "start",
                       drop = TRUE]

  names(opt_params) <- opt_param_names

  if (fit_decision %in% "continue") {

    lower_params <- par_DF[which(par_DF$optimize_param), "lower_bound",
                           drop = TRUE]
    names(lower_params) <- opt_param_names

    upper_params <- par_DF[which(par_DF$optimize_param), "upper_bound",
                           drop = TRUE]
    names(upper_params) <- opt_param_names




    # Now call optimx::optimx() and do the fit
    optimx_out <- suppressWarnings(
      tryCatch(
        expr = {
          tmp <- do.call(
            optimx::opm,
            args = c(
              list(par = opt_params,
                   fn = log_likelihood,
                   lower = lower_params,
                   upper = upper_params),
              # method and control
              settings_optimx,
              # ... additional args to log_likelihood
              list(
                const_params = const_params,
                data = data,
                data_sigma_group = data$data_sigma_group,
                modelfun = modelfun,
                dose_norm = dose_norm,
                log10_trans = log10_trans,
                negative = TRUE,
                force_finite = TRUE,
                suppress.messages = suppress.messages
              ) # end list()
            ) # end args = c()
          ) # end do.call
          tmp$method <- rownames(tmp)
          tmp <- merge(tmp, as.data.frame(attr(tmp, "details")))
          tmp
        },
        error = function(err) {
          method <- rlang::eval_tidy(settings_optimx$method)
          tmp <- c(rep(NA_real_, length(opt_params)),
                   rep(NA_real_, 4),
                   -9999,
                   rep(NA, 2),
                   NA_real_)

          names(tmp) <- c(names(opt_params),
                          "value", "fevals", "gevals", "hevals",
                          "convergence",
                          "kkt1", "kkt2",
                          "xtime")

          tmp <- data.frame(as.list(tmp)) |>
            dplyr::slice(rep(seq_len(dplyr::n()),
                             each = length(method)))

          rownames(tmp) <- method
          tmp$method <- rownames(tmp)

          details <- as.data.frame(
            cbind(method = as.list(method),
                  ngatend = as.list(rep(NA_real_, length(method))),
                  nhatend = as.list(rep(NA_real_, length(method))),
                  hev = as.list(rep(NA_real_, length(method))),
                  message = as.list(rep(err$message,
                                        length(method)))
            )
          )
          rownames(details) <- method
          tmp <- merge(tmp, details)

          return(tmp)
        }) # end tryCatch()
    ) # end suppressWarnings()
    out <- optimx_out
  } else { # if status for this model was "abort", then abort fit and return NULL
    method <- rlang::eval_tidy(settings_optimx$method)
    tmp <- c(rep(NA_real_, length(opt_params)),
             rep(NA_real_, 4),
             9999,
             rep(NA, 2),
             NA_real_)
    names(tmp) <- c(names(opt_params),
                    "value", "fevals", "gevals", "hevals",
                    "convergence",
                    "kkt1", "kkt2",
                    "xtime")
    tmp <- data.frame(as.list(tmp)) |>
      dplyr::slice(rep(seq_len(dplyr::n()),
                       each = length(method)))
    rownames(tmp) <- method
    tmp$method <- rownames(tmp)

    attr(tmp, "maximize") <- FALSE
    attr(tmp, "npar") <- length(opt_params)
    attr(tmp, "follow.on") <- FALSE
    details <- as.data.frame(
      cbind(method = as.list(method),
            ngatend = as.list(rep(NA_real_, length(method))),
            nhatend = as.list(rep(NA_real_, length(method))),
            hev = as.list(rep(NA_real_, length(method))),
            message = as.list(rep("Fit status was 'abort'",
                                  length(method)))
      )
    )
    rownames(details) <- method
    tmp <- merge(tmp, details)

    out <- tmp
  }

  return(out)
}
