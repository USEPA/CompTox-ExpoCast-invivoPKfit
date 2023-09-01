#' Fit a single group of data
#'
#' @param data A single group of data
#' @param sigma_DF sigma_DF for a single group of data
#' @param par_DF par_DF for a single group of data
#' @param this_model Name of the `pk_model` object to fit
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
                      suppress.messages){
  #Rowbind par_DF and sigma_DF
  par_DF <- dplyr::bind_rows(par_DF,
                             sigma_DF)

  #get params to be optimized with their starting points
  opt_params <- par_DF %>%
    dplyr::filter(optimize_param %in% TRUE) %>%
    dplyr::pull(start)
  names(opt_params) <- par_DF %>%
    dplyr::filter(optimize_param %in% TRUE) %>%
    dplyr::pull(param_name)

  if (fit_decision %in% "continue") {

    #param lower bounds (only on params to be optimized)
    lower_params <- par_DF %>%
      dplyr::filter(optimize_param %in% TRUE) %>%
      dplyr::pull(lower_bound)
    names(lower_params) <- par_DF %>%
      dplyr::filter(optimize_param %in% TRUE) %>%
      dplyr::pull(param_name)

    #param upper bounds (only on params to be optimized)
    upper_params <- par_DF %>%
      dplyr::filter(optimize_param %in% TRUE) %>%
      dplyr::pull(upper_bound)
    names(upper_params) <- par_DF %>%
      dplyr::filter(optimize_param %in% TRUE) %>%
      dplyr::pull(param_name)

    #get params to be held constant, if any
    if(any(par_DF$optimize_param %in% FALSE &
           par_DF$use_param %in% TRUE)){
      const_params <- par_DF %>%
        dplyr::filter(optimize_param %in% FALSE &
                        use_param %in% TRUE) %>%
        dplyr::pull(start)

      names(const_params) <- par_DF %>%
        dplyr::filter(optimize_param %in% FALSE &
                        use_param %in% TRUE) %>%
        dplyr::pull(param_name)
    }else{
      const_params <- NULL
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
            settings_optimx,
            #... additional args to log_likelihood
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
            ) #end list()
          ) #end args = c()
        ) #end do.call
        tmp$method <- rownames(tmp)
        tmp
      },
      error = function(err){
        method <- settings_optimx$method
        tmp <- data.frame(c(rep(NA_real_,
                              length(opt_params)),
                          rep(NA_real_, 4),
                          -9999,
                          rep(NA, 2),
                          NA_real_))

        names(tmp) <- c(names(opt_params),
                        "value", "fevals", "gevals", "niter",
                        "convcode",
                        "kkt1", "kkt2",
                        "xtime")
        rownames(tmp) <- method
        tmp$method <- rownames(tmp)

        attr(tmp, "maximize") <- FALSE
        attr(tmp, "npar") <- length(opt_params)
        attr(tmp, "follow.on") <- FALSE

        details <- cbind(method = as.list(method),
                                           ngatend = as.list(rep(NA_real, length(method))),
                                           nhatend = as.list(rep(NA_real, length(method))),
                                           hev = as.list(rep(NA_real, length(method))),
                                           message = as.list(rep(err$message,
                                                                 length(method)
                                                                 )
                                                             )
                                           )
        rownames(details) <- method
        attr(tmp, "details") <- details

        class(tmp) <- c("optimx",
                        class(tmp))

        return(tmp)
      }) #end tryCatch()
    ) #end suppressWarnings()

    # Need to convert units back

    out <- optimx_out

  }else{ #if status for this model was "abort", then abort fit and return NULL
    method <- rlang::eval_tidy(settings_optimx$method)
    tmp <- c(rep(NA_real_,
                            length(opt_params)),
                        rep(NA_real_, 4),
                        -9999,
                        rep(NA, 2),
                        NA_real_)
    names(tmp) <- c(names(opt_params),
                    "value", "fevals", "gevals", "niter",
                    "convcode",
                    "kkt1", "kkt2",
                    "xtime")
    tmp <- data.frame(as.list(tmp)) %>%
      dplyr::slice(rep(1:dplyr::n(),
                       each = length(method)))
    rownames(tmp) <- method
    tmp$method <- rownames(tmp)

    attr(tmp, "maximize") <- FALSE
    attr(tmp, "npar") <- length(opt_params)
    attr(tmp, "follow.on") <- FALSE
    # method <- settings_optimx$method
    details <- cbind(method = as.list(method),
                     ngatend = as.list(rep(NA_real_, length(method))),
                     nhatend = as.list(rep(NA_real_, length(method))),
                     hev = as.list(rep(NA_real_, length(method))),
                     message = as.list(rep("Fit status was 'abort'",
                                           length(method)
                     )
                     )
    )
    rownames(details) <- method
    attr(tmp, "details") <- details

    class(tmp) <- c("optimx",
                    class(tmp))

    out <- tmp
  }

  return(out)
}
