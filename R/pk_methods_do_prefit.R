#' Do pre-fitting
#'
#' Do pre-fit calculations and checks
#'
#' This function does the following:
#' \itemize{
#' \item Based on the error group in `pk_groups` and the pre-processed data, determines the number of residual standard deviations ("sigmas") hyperparameters to be estimated.
#' \item Determines which "sigma" hyperparameter corresponds to each observation in the data.
#' \item Calculates lower/upper bounds and starting guesses for each "sigma" hyperparameter
#' \item For each model in `stat_model`, calls its `params_fun`, the function that, based on the data, determines whether to optimize each model parameter, and calculates lower/upper bounds and starting guesses for each model parameter to be optimized. Only non-excluded observations are passed to each model's `params_fun`.
#' }
#'
#' Lower bounds for each "sigma" hyperparameter are set to `sqrt(.Machine$double_eps)`.
#'
#' Upper bounds for each "sigma" hyperparameter are calculated as the standard
#' deviation of observations in the corresponding error SD group (see
#' [combined_sd()]), with any specified transformations applied
#' (dose-normalization and/or log10-transformation). If the combined SD is
#' non-finite or less than the sigma lower bound, then the maximum concentration
#' is used as an upper bound; if this still returns a non-finite value or a
#' value less than the lower bound, then a constant value of 1000 is
#' substituted.
#'
#' The starting guess for each "sigma" hyperparameter is one-tenth of the upper bound.
#'
#' If there are less detected observations than timepoints, or if there are
#' parameters necessary for model fitting that have missing values,
#' these models will not be fit.
#'
#'
#' @inheritParams do_preprocess.pk
#' @return The same `pk` object, but with a new element `prefit`, containing the
#'   results of pre-fit calculations and checks for each model and for the error
#'   model.
#' @export
#' @author Caroline Ring
do_prefit.pk <- function(obj, ...) {

    objname <- deparse(substitute(obj))
    status <- obj$status
    if (status >= status_prefit) {
      warning(objname,
              " current status is ",
    status,
    ". do_prefit.pk() will reset its status to ",
    status_prefit,
    ". Any results from later workflow stages will be lost."
      )
    }
    cli::cli_par()

    # if preprocessing not already done, do it
    if (obj$status < status_preprocess) {
      obj <- do_preprocess(obj)
    }

    if (obj$status < status_data_info) {
      obj <- do_data_info(obj)
    }

    suppress_messages <- obj$pk_settings$preprocess$suppress.messages

    data <- get_data.pk(obj)
    data_grp <- get_data_group.pk(obj)
    data_grp_vars <- get_data_group.pk(obj, as_character = TRUE)
    data_group_keys <- data |> dplyr::distinct(DATA_GROUP_ID, !!!data_grp)

    if (suppress_messages %in% FALSE) {
      cli_inform("do_prefit.pk(): Assigning error SD groups to all observations")
    }

    error_group <- get_error_group.pk(obj)
    # get the error model/group from `pk_groups`, which defines the number of sigmas that will need to be optimized
    # count the number of unique combinations of vars in obj$stat_error_model$error_group
    unique_groups <- unique(rlang::eval_tidy(expr = error_group, data = data))
    # Assign a factor variable denoting sigma group to each observation in obj$data
    # This tells us which sigma applies to which observation
    data_sigma_group <- interaction(
      lapply(error_group, rlang::eval_tidy, data = data),
    sep = "_"
    )

    # but set to NA for any excluded observations, and drop any levels that occur only for excluded observations
    data_sigma_group[data$exclude %in% TRUE] <- NA_character_
    data_sigma_group <- droplevels(data_sigma_group)

    # add data_sigma_group to data
    obj$data$data_sigma_group <- data_sigma_group
    data <- get_data.pk(obj)

    # get bounds and starting points for each error sigma to be fitted
    if (suppress_messages %in% FALSE) {
      cli_inform(c(
        paste("do_prefit.pk(): ",
              "Getting bounds and starting guesses for each error SD to be fitted"
        )
      ))
    }

    # Set a value to square root of lowest possible value 'x' where 1+x != 1
    sigma_lower <- sqrt(.Machine$double.eps)

    # Add the data_sigma_group and filter out all excluded values
    sigma_DF <- data |>
      dplyr::mutate(data_sigma_group = data_sigma_group) |>
      dplyr::filter(exclude %in% FALSE) |>
      # temporarily undo log10-trans, if it has been used
      # this is because combined_sd() requires NON log transformed concs
      # but we want to keep dose-normalization if it has been applied
      # so we un-log Conc_trans
      dplyr::mutate(
        Conc_tmp = dplyr::if_else(rep(obj$scales$conc$log10_trans %in% TRUE,
                                      NROW(Conc_trans)),
                                  10^Conc_trans,
                                  Conc_trans),
        Conc_SD_tmp = dplyr::if_else(rep(obj$scales$conc$log10_trans %in% TRUE,
                                         NROW(Conc_trans)),
                                     10^Conc_SD_trans,
                                     Conc_SD_trans),
        Conc_tmp.Units = dplyr::if_else(rep(obj$scales$conc$log10_trans %in% TRUE,
                                            NROW(Conc_trans)),
                                        gsub(pattern = "\\)$",
                                             replacement = "",
                                             x = gsub(pattern = "^log10\\(",
                                                      replacement = "",
                                                      x = Conc_trans.Units)),
                                        Conc_trans.Units)
      )

    # Set values for sigma upper/lower-bounds and start
    sigma_DF <- sigma_DF |>
      dplyr::group_by(!!!error_group) |>
      dplyr::summarise(param_name = paste("sigma",
                                          unique(data_sigma_group),
                                          sep = "_"),
                       param_units = unique(Conc_tmp.Units),
                       optimize_param = TRUE,
                       use_param = TRUE,
                       lower_bound = sigma_lower,
                       upper_bound = combined_sd(
                         group_mean = Conc_tmp,
                         group_sd = Conc_SD_tmp,
                         group_n = N_Subjects,
                         unbiased = TRUE,
                         na.rm = TRUE,
                         log10 = obj$scales$conc$log10_trans)) |>
      dplyr::mutate(
        upper_bound = dplyr::if_else(
          !is.finite(upper_bound) | upper_bound <= lower_bound,
          1000,
          upper_bound),
        start = 0.1 * upper_bound) |>
      as.data.frame()

    # assign rownames to sigma_DF
    rownames(sigma_DF) <- sigma_DF$param_name

    # assign sigma_DF to the `pk` object
    obj$prefit$sigma_DF <- sigma_DF

    if (suppress_messages %in% FALSE) {
      cli_inform(c(
        paste("do_prefit.pk(): Getting bounds and starting guesses",
              "for all model parameters to be fitted"
        )
      ))
    }

    # for each model to be fitted:
    par_DF_out <- sapply(
      names(obj$stat_model),
      function(this_model) {
        # get parameters to be optimized, bounds, and starting points
        # by evaluating params_fun for this stat_model
        # pass it only the non-excluded observations
        # NOTE: data is passed on so data parameters can be set within the function
        base::split.data.frame(data, f = data[["DATA_GROUP_ID"]]) |>
          lapply(
            function(x) {
              fx <- subset(x, exclude %in% FALSE)
              fx <- do.call(
                what = obj$stat_model[[this_model]]$params_fun,
                args = append(list(data = fx), obj$stat_model[[this_model]]$params_fun_args)
              )
              # Output must have the DATA_GROUP_ID
              # cbind(DATA_GROUP_ID = unique(x[["DATA_GROUP_ID"]]), fx) |>
              #   tibble::remove_rownames()
              tibble::remove_rownames(fx)
            }
          ) |>
          dplyr::bind_rows(.id = "DATA_GROUP_ID") |>
          dplyr::mutate(DATA_GROUP_ID = as.integer(DATA_GROUP_ID))
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )

    par_DF_out <- dplyr::bind_rows(par_DF_out, .id = "model") |>
      dplyr::left_join(
        data_group_keys,
        by = "DATA_GROUP_ID",
        relationship = "many-to-one"
      ) |>
      dplyr::relocate(!!!data_grp, .after = "DATA_GROUP_ID")

    fit_check_out <- sapply(
      names(obj$stat_model),
      function(this_model) {
        # check whether there are enough observations to optimize the requested parameters plus sigmas
        # number of parameters to optimize
        if (suppress_messages %in% FALSE) {
          cli_inform(
            "do_prefit.pk(): Checking whether sufficient observations to fit models."
          )
        }
        # Are any parameters used initialized to NA?
        n_par_DF <- par_DF_out |>
          dplyr::filter(model %in% this_model) |>
          dplyr::group_by(!!!data_grp) |>
          dplyr::summarise(
            n_par = sum(optimize_param),
            used_par_na = any(use_param & is.na(start)
            )
          )

        n_sigma_DF <- sigma_DF |>
          dplyr::group_by(!!!data_grp) |>
          dplyr::summarise(n_sigma = sum(optimize_param))


        n_detect_DF <- get_data_summary(obj) |>
          dplyr::group_by(!!!data_grp) |>
          dplyr::summarise(n_detect = sum(n_detect))

        # merge all of these together
        fit_check_DF <- dplyr::inner_join(
          dplyr::inner_join(n_par_DF,
                            n_sigma_DF,
                            by = data_grp_vars),
          n_detect_DF,
          by = data_grp_vars
        )

        # get fit decision & reasoning
        fit_check_DF <- fit_check_DF |>
          dplyr::mutate(
            n_par_opt = n_par + n_sigma,
            fit_decision = ifelse(n_par_opt < n_detect & !used_par_na,
                                  "continue",
                                  "abort"),
            fit_reason = dplyr::case_when(
              used_par_na ~ "Some parameters necessary for model fitting are NA.",
              n_par_opt < n_detect ~ paste(
                "Number of parameters to estimate is ",
                "less than number of non-excluded detected observations"),
              .default = paste(
                "Number of parameters to estimate is ",
                "greater than or equal to number of non-excluded detected observations")
            )
          ) |> as.data.frame()


        fit_check_DF
      },
      simplify = FALSE,
      USE.NAMES = TRUE)

    fit_check_out <- dplyr::bind_rows(fit_check_out, .id = "model")

    par_DF_out <- par_DF_out |>
      dplyr::mutate(
        lower_bound = ifelse(
          is.na(lower_bound) & !is.na(start) & optimize_param,
          start, lower_bound
        ),
        upper_bound = ifelse(
          is.na(upper_bound) & !is.na(start) & optimize_param,
          start, upper_bound
        )
      )

    obj$prefit$par_DF <- par_DF_out
    obj$prefit$fit_check <- fit_check_out

    obj$status <- status_prefit # prefit complete

    return(obj)
}
