#' Hyperparameter fitting
#'
#' Fit hyperparameter sigma for a `pk` object with pre-calculated model predictions
#'
#' This function estimates the hyperparameter \eqn{\sigma} from
#' a data.frame of pre-calculated model predictions,
#' using numerical optimization implemented in [optimx::optimx()]. The
#' optimization is done by maximizing the log-likelihood function implemented in
#' [log_likelihood()]. Only the non-excluded observations are used.
#'
#'
#' @param obj A [pk] object.
#' @param preds A data.frame similar to the results from [predict.pk()] which
#' contains pre-calculated predictions in addition to the concentration over time
#' values which can be obtained from [get_data.pk()].
#' @param pred_col A character vector with the name of the column with predictions.
#' @param k Default 2. The `k` parameter in the log-likelihood formula (see
#' Details). Must be named if used.
#' @param ... Additional arguments. Not currently in use.
#' @return The same [pk] object, with new element beginning with `ext_fit`
#' containing a list of the prediction data used as input and two data.frames with
#' optimized sigma values and AICs for those predictions per `data_group`.
#' @export
#' @author Gilberto Padilla Mercado
#'

fit_sigma.pk <- function(obj, preds, pred_col, k = 2, ...) {
  # check status
  objname <- deparse(substitute(obj))
  status <- obj$status

  if (status < status_prefit) {
    stop(
      objname, " current status is ",
      status, ". Needs to be at least prefit. ",
      "(status >= 3)."
    )
  }
  # Get data groups for referencing
  data <- get_data(obj)

  data_grp <- get_data_group.pk(obj)
  data_grp_vars <-  get_data_group.pk(obj, as_character = TRUE)

  # Get the scaling transformations
  this_settings_optimx <- get_settings_optimx(obj)
  dose_norm <- obj$scales$conc$dose_norm
  log10_trans <- obj$scales$conc$log10_trans

  # Needs more flexible naming checks, fix in future commits
  req_names <- c(
    "Chemical", "Species", "Dose", "Route", "Media", "Reference",
    "N_Subjects", "pLOQ", "Time", "Conc", "Conc_trans", "Conc_SD",
    "data_sigma_group", "Detect", "exclude"
  )

  data_names <- names(preds)
  if (!all(req_names %in% data_names)) {
    cli::cli_abort(
      "Please make sure that there are the following columns names: {paste(req_names)}"
    )
  }

  # Must be the name of a single column (for now)
  stopifnot(length(pred_col) == 1)

  if (!(pred_col %in% names(preds))) {
    cli::cli_abort("Prediction column not in prediction data.frame.")
  }


  # Inner join the predictions and observations (there should be a lot of overlap)
  sigma_df <- obj$prefit$stat_error_model$sigma_DF |>
    dplyr::select(!!!data_grp, "data_sigma_group" = param_name,
                  lower_bound, upper_bound, start) |>
    dplyr::distinct()

  data_preds <- dplyr::inner_join(data, preds) |>
    dplyr::select(c(req_names, pred_col)) |>
    dplyr::mutate(data_sigma_group = paste0("sigma_", data_sigma_group)) |>
    dplyr::inner_join(sigma_df) |>
    dplyr::rename("Conc_est" = !!pred_col)

  # Split data
  sigma_group_fct <- factor(data_preds$data_sigma_group)
  sp_data <- split(data_preds, sigma_group_fct)
  sp_data # per split needs to have opt_params, lower_bounds, upper_bounds

  setup_df <- lapply(sp_data, \(x) {
    param_start <-   unique(x$start)
    param_lower <-   unique(x$lower_bound)
    param_upper <-   unique(x$upper_bound)

    names(param_start) <- unique(x$data_sigma_group)
    names(param_lower) <- unique(x$data_sigma_group)
    names(param_upper) <- unique(x$data_sigma_group)

    list(data = x,
         param_start = param_start,
         param_lower = param_lower,
         param_upper = param_upper
    )
    }
    )

  optimized_df <- lapply(
    setup_df,
    \(x) {
      tryCatch(
        expr = {
          tmp <- do.call(
            optimx::optimx,
            args = c(
              list(par = x$param_start,
                   fn = log_likelihood,
                   lower = x$param_lower,
                   upper = x$param_upper
              ),
              # method and control
              this_settings_optimx,
              # ... additional args to log_likelihood
              list(
                const_params = NULL,
                data = x$data,
                data_sigma_group = x$data$data_sigma_group,
                modelfun = pred_col,
                dose_norm = dose_norm,
                log10_trans = log10_trans,
                negative = TRUE,
                force_finite = TRUE,
                includes_preds = TRUE,
                suppress.messages = TRUE
              ) # end list()
            ) # end args = c()
          ) # end do.call
          tmp$method <- rownames(tmp)
          tmp <- merge(tmp, as.data.frame(attr(tmp, "details")))
          tmp
        },
        error = function(err) {
          method <- rlang::eval_tidy(this_settings_optimx$method)
          tmp <- c(rep(NA_real_, length(x$param_start)),
                   rep(NA_real_, 4),
                   -9999,
                   rep(NA, 2),
                   NA_real_)

          names(tmp) <- c(names(x$param_start),
                          "value", "fevals", "gevals", "niter",
                          "convergence",
                          "kkt1", "kkt2",
                          "xtime")

          tmp <- data.frame(as.list(tmp)) |>
            dplyr::slice(rep(seq_len(dplyr::n()), each = length(method)))

          rownames(tmp) <- method
          tmp$method <- rownames(tmp)

          details <- as.data.frame(
            cbind(method = as.list(method),
                  ngatend = as.list(rep_len(NA_real_, length(method))),
                  nhatend = as.list(rep_len(NA_real_, length(method))),
                  hev = as.list(rep_len(NA_real_, length(method))),
                  message = as.list(rep_len(err$message, length(method)))
            )
          )
          rownames(details) <- method
          tmp <- merge(tmp, details)

          return(tmp)
        })
  })

  optimized_df <- lapply(
    optimized_df,
    \(x) {
      x |>
        tidyr::pivot_longer(cols = dplyr::starts_with("sigma_"),
                            names_to = "hyperparam_name",
                            values_to = "hyperparam_value")
    }
  )

  final_df <- right_join(
    sigma_df,
    dplyr::bind_rows(optimized_df) |>
      dplyr::mutate(value = -1 * value, # Negative values used in minimization
             ngatend = as.numeric(ngatend),
             nhatend = as.numeric(nhatend),
             hev = as.numeric(hev),
             model = pred_col
             ) |>
      dplyr::relocate(dplyr::starts_with("hyperparam_")),
    by = dplyr::join_by("data_sigma_group" == "hyperparam_name"))

  AIC_df <- final_df |>
    dplyr::select(!!!data_grp,
                  model, method,
                  data_sigma_group,
                  value) |>
    dplyr::group_by(!!!data_grp, model, method) |>
    dplyr::summarize(ll = sum(value),
                     npar = n(),
                     data_sigma_group = c(data_sigma_group)
    ) |>
    dplyr::group_by(!!!data_grp, model, method) |>
    dplyr::mutate(AIC = (k * .data$npar) - (2 * .data$ll))
  # Usually only one parameter is optimized, but sometimes there is more than one
  # sigma value per data group.

  obj[[paste0("ext_fits.", pred_col)]] <- list(
    sigma_fits = final_df,
    sigma_AIC = AIC_df,
    input_preds = preds
  )


  return(obj)
}
