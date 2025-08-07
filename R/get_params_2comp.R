#' Get parameters for 2-compartment model
#'
#' Get parameters for 2-compartment model and determine whether each is to be
#' estimated from the data
#'
#' The full set of model parameters for the 2-compartment model includes `V1`,
#' `kelim`, `k12`, `k21`, `kgutabs`,`Fgutabs`, and `Rblood2plasma`. Whether each one can be
#' estimated from the data depends on what routes of administration are included
#' in the data.
#'
#' @section IV data, no oral data:
#'
#' If IV dosing data are available, but no oral dosing data are available, then
#' only the parameters `V1`, `kelim`, `k12`, and `k21` will be estimated from
#' the data. The parameters `kgutabs` and `Fgutabs` cannot be estimated from IV
#' data alone.
#'
#' @section Oral data, no IV data:
#'
#' If oral dosing data are available, but no IV dosing data are available, then
#' the parameters `kelim`, `k12`, `k21`, and `kgutabs` will be estimated from
#' the data. However, the parameters `Fgutabs` and `V1` cannot be identified
#' separately. From oral data alone, only the ratio `Fgutabs/V1` can be
#' identified. This ratio is represented by a single parameter named
#' `Fgutabs_V1`. `Fgutabs` and `V1` will not be optimized, but `Fgutabs_V1` will
#' be optimized, along with `kelim`, `k12`, `k21`, and `kgutabs`.
#'
#' @section Oral data and IV data:
#' If both oral and IV dosing data are available, then `V1`, `kelim`, `k12`,
#' `k21`, `kgutabs`, and `Fgutabs` will all be estimated from the data.
#'
#' @inheritSection get_params_flat Blood and plasma data
#' @inheritSection get_params_flat Only one of blood or plasma data
#'
#' @section Default lower and upper bounds for each parameter:
#' \subsection{Default lower and upper bounds for time constants `kelim`, `kgutabs`,
#' `k12`, and `k21`.}{
#' Default bounds for time constants `kelim` and `kgutabs` are set based on
#' the time scale of the available data.
#'
#' The lower bounds are based on the assumption that elimination, absorption,
#' and distribution are very slow compared to the time scale of the study.
#' Specifically, the lower bounds assume thatelimination, absorption, and
#' distribution half-lives are twice as long as the duration of the available
#' study data, or `2*max(Time_trans)`. Under this assumption, the corresponding
#' elimination, absorption, and distribution time constants would be
#' `log(2)/(2*max(Time_trans))`. Therefore, the default lower bounds for
#' `kelim`, `kgutabs`, `k12`, and `k21` are `log(2)/(2*max(Time_trans))`.
#'
#' Upper bounds are based on the opposite assumption: that elimination,
#' absorption, and distribution are very fast compared to the time scale of the
#' study. Specifically, the upper bounds assume that the elimination,
#' absorption, and distribution half-lives are half as long as the time of the
#' first observation after time 0, or `0.5*min(Time_trans[Time_trans>0])`. Under
#' this assumption, the correspondingelimination, absorption, and distribution
#' time constants would be `log(2)/(0.5*min(Time_trans[Time_trans>0]))`.
#' Therefore, the default lower bounds for `kelim`, `kgutabs`, `k12`, and `k21`
#' are `log(2)/(0.5*min(Time_trans[Time_trans>0]))`.
#' }
#' \subsection{Default lower and upper bounds for `V1`}{
#' By default, the lower bound for `V1` is 0.01, and the upper bound for
#' `V1` is 100. These values were chosen based on professional judgment.
#' }
#' \subsection{Default lower and upper bounds for `Fgutabs`}{
#' By default, the lower bound for `Fgutabs` is 0.0, and the upper bound for
#' `Fgutabs` is 1. These are simply the bounds of the physically-meaningful
#' range for a fraction.
#' }
#' \subsection{Default lower and upper bounds for `Fgutabs_V1`}{
#' By default, the lower bound for the ratio `Fgutabs_V1` is 0.01, and the
#' upper bound is 100. These values were chosen based on professional judgment.
#' }
#' \subsection{Default lower and upper bounds for `Rblood2plasma`}{
#' By default, the lower bound for the blood:plasma partition coefficient
#' `Rblood2plasma` is 0.01, and the upper bound is 100. These values were chosen
#' based on professional judgment.
#' }
#'
#'
#' @section Starting values for each parameter:
#'
#' Starting values for each parameter (starting guesses for the numerical
#' optimizer) are derived from the data using [get_starts_2comp()].
#'
#' If the starting values returned by [get_starts_2comp()] fall outside the
#' bounds for any parameter(s), then the starting value will be reset to a value
#' halfway between the lower and upper bounds for that parameter.
#'
#' @inheritParams get_params_flat
#'
#' @return A `data.frame`with the following variables:
#' \itemize{
#' \item `param_name`: Character: Names of the model parameters
#' \item `param_units`: Character: Units of the model parameters
#' \item `optimize_param`: TRUE if each parameter is to be estimated from the data; FALSE otherwise
#' \item `use_param`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
#' \item`lower_bounds`: Numeric: The lower bounds for each parameter
#' \item `upper_bounds`: Numeric: The upper bounds for each parameter
#' \item `start`: Numeric: The starting guesses for each parameter
#' }
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family 2-compartment model functions
#' @family get_params functions
#' @family built-in model functions

get_params_2comp <- function(
    data,
    lower_bound = NULL,
    upper_bound = NULL,
    param_units = alist(
      kelim = paste0("1/", unique(Time_trans.Units)),
      V1 = paste0("(", unique(Dose.Units), ")/(", unique(Conc.Units), ")"),
      k21 = paste0("1/", unique(Time_trans.Units)),
      k12 = paste0("1/", unique(Time_trans.Units)),
      Fgutabs = "unitless fraction",
      kgutabs = paste0("1/", unique(Time_trans.Units)),
      Fgutabs_V1 = paste0("(", unique(Conc.Units), ")/(", unique(Dose.Units), ")"),
      Rblood2plasma = "unitless ratio"
    ),
    ...) {
  # param names
  param_name <- c(
    "kelim",
    "V1",
    "k21",
    "k12",
    "Fgutabs",
    "kgutabs",
    "Fgutabs_V1",
    "Rblood2plasma"
  )
  # Collect extra parameters
  dots <- list(...)

  # Default lower bounds, to be used in case the user specified a non-default
  # value for the `lower_bound` argument, but did not specify expressions for all
  # parameters. Any parameters not specified in the `lower_bound` argument will
  # take their default lower bounds defined here. This should be the same as the
  # default value for the `lower_bound` argument.
  lower_bound_default <- alist(
    kelim = log(2) / (2 * max(Time_trans)),
    k12 = log(2) / (2 * max(Time_trans)),
    k21 = log(2) / (2 * max(Time_trans)),
    V1 = 0.01,
    Fgutabs = 1E-3,
    kgutabs = log(2) / (2 * max(Time_trans)),
    Fgutabs_V1 = 1E-5,
    Rblood2plasma = 1e-2
  )

  # which parameters did not have lower bounds specified in the `lower_bound`
  # argument?
  lower_bound_missing <- base::setdiff(
    names(lower_bound_default),
    names(lower_bound)
  )
  # fill in the default lower bounds for any parameters that don't have them
  # defined in the `lower_bound` argument
  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]

  # Default upper bounds, to be used in case the user specified a non-default
  # value for the `upper_bound` argument, but did not specify expressions for all
  # parameters. Any parameters not specified in the `upper_bound` argument will
  # take their default upper bounds defined here. This should be the same as the
  # default value for the `upper_bound` argument.
  upper_bound_default <- alist(
    kelim = log(2) / (0.5 * min(Time_trans[Time_trans > 0])),
    k12 = log(2) / (0.5 * min(Time_trans[Time_trans > 0])),
    k21 = log(2) / (0.5 * min(Time_trans[Time_trans > 0])),
    V1 = 100.0,
    Fgutabs = 1.0,
    kgutabs = log(2) / (0.5 * min(Time_trans[Time_trans > 0])),
    Fgutabs_V1 = 1e2,
    Rblood2plasma = 100
  )

  # which parameters did not have upper bounds specified in the `upper_bound`
  # argument?
  upper_bound_missing <- base::setdiff(
    names(upper_bound_default),
    names(upper_bound)
  )
  # fill in the default upper bounds for any parameters that don't have them
  # defined in the `upper_bound` argument
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]

  # initialize optimization: start with optimize = TRUE for all params
  optimize_param <- rep_len(TRUE, length(param_name))

  # initialize whether each param is used: start with use = TRUE for all params
  use_param <- rep_len(TRUE, length(param_name))

  has_oral <- "oral" %in% data$Route
  has_iv <- "iv" %in% data$Route

  if ("pars_to_optimize" %in% names(dots)) {
    pars_to_optimize <- dots[["pars_to_optimize"]]
    pars_to_hold <- setdiff(param_name, pars_to_optimize)
    optimize_param[param_name %in% pars_to_hold] <- FALSE

    if (xor(has_oral, has_iv) &&
      "Fgutabs" %in% pars_to_optimize &&
      !"Fgutabs_V1" %in% pars_to_optimize) {
      cli::cli_warn("`Fgutabs` will not be optimized due to contraints in data: Route")
    }
    if (length(pars_to_hold) == length(param_name)) {
      cli::cli_warn("No parameters to optimize specified.")
    }
  }

  # Three mutually exclusive conditions based on Route:
  if (has_oral && has_iv) {
    # if both oral and IV data, Fgutabs and V1 can be fit separately, so turn off Fgutabs_V1
    optimize_param[param_name %in% "Fgutabs_V1"] <- FALSE
    use_param[param_name %in% "Fgutabs_V1"] <- FALSE
  }

  if (has_iv && !has_oral) {
    # if no oral data, can't fit kgutabs, Fgutabs, or Fgutabs_V1,
    # and they won't be used.
    optimize_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_V1")] <- FALSE
    use_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_V1")] <- FALSE
  }

  if (has_oral && !has_iv) {
    # if oral ONLY:
    # cannot fit Fgutabs and V1 separately, so turn them off.
    optimize_param[param_name %in% c("Fgutabs", "V1")] <- FALSE
    use_param[param_name %in% c("Fgutabs", "V1")] <- FALSE
  }


  # if both "blood" and "plasma" are not in data,
  # then Rblood2plasma will not be optimized
  if (!(all(c("blood", "plasma") %in% data$Media))) {
    optimize_param[param_name %in% "Rblood2plasma"] <- FALSE
  }

  # get param units based on data
  param_units_vect <- sapply(param_units,
    rlang::eval_tidy,
    data = data,
    simplify = TRUE,
    USE.NAMES = TRUE
  )
  param_units_vect <- param_units_vect[param_name]


  lower_bound_vect <- sapply(lower_bound,
    rlang::eval_tidy,
    data = data,
    simplify = TRUE,
    USE.NAMES = TRUE
  )
  lower_bound_vect <- lower_bound_vect[param_name]

  upper_bound_vect <- sapply(upper_bound,
    rlang::eval_tidy,
    data = data,
    simplify = TRUE,
    USE.NAMES = TRUE
  )
  upper_bound_vect <- upper_bound_vect[param_name]

  par_DF <- data.frame(
    "param_name" = param_name,
    "param_units" = param_units_vect,
    "optimize_param" = optimize_param,
    "use_param" = use_param,
    "lower_bound" = lower_bound_vect,
    "upper_bound" = upper_bound_vect
  )

  par_DF <- get_starts_2comp(data = data, par_DF = par_DF, ...)

  # check to ensure starting values are within bounds
  # if not, then replace them by a value halfway between bounds
  start_low <- (par_DF$start < par_DF$lower_bound) %in% TRUE
  start_high <- (par_DF$start > par_DF$upper_bound) %in% TRUE
  start_nonfin <- !is.finite(par_DF$start)

  par_DF[
    start_low | start_high | start_nonfin,
    "start"
  ] <- rowMeans(
    cbind(par_DF[
      start_low | start_high | start_nonfin,
      c(
        "lower_bound",
        "upper_bound"
      )
    ])
  )

  par_DF$start <- as.list(par_DF$start)

  return(par_DF)
}
