#' Get parameters for 1-compartment model with excretion
#'
#' Get parameters for 1-compartment model and determine whether each is to be
#' estimated from the data. This is intended for use with collected radiolabeling
#' data.
#'
#' The full set of model parameters for the 1-compartment model includes
#' `Clint`,  `Fup`, `kgutabs`, `Fgutabs`, `Vdist`,
#' `Q_totli`, `Q_gfr`, `Q_alv`, `Kblood2air`, and `Rblood2plasma`.
#' Since these are directly taken from `httk`, many of the usual route-dependent
#' estimations of these parameters are not necessary.
#'
#' @inheritSection get_params_1comp_cl Constant parameters from `httk`
#' @inheritSection get_params_1comp_cl IV data, no oral data
#' @inheritSection get_params_1comp_cl Oral data, no IV data
#' @inheritSection get_params_1comp_cl Oral data and IV data
#' @inheritSection get_params_flat Blood and plasma data
#' @inheritSection get_params_flat Only one of blood or plasma data
#'
#' @section Fraction recovered `Frec`:
#'
#' This is a parameter that is estimated by taking the maximum value
#' for the excreta medium measurements (which as a reminder is a cumulative
#' measurement) and dividing this value by the dose (which should be 100% of
#' what goes into the body). For oral doses this value is divided by `Dose * Fgutabs`.
#'
#' @section Default lower and upper bounds for each parameter:
#' \subsection{Default lower and upper bounds for `kelim` and `kgutabs`}{
#' Default bounds for time constants `kelim` and `kgutabs` are set based on
#' the time scale of the available data.
#'
#' The lower bounds are based on the assumption that elimination and absorption
#' are very slow compared to the time scale of the study. Specifically, the
#' lower bounds assume that elimination and absorption half-lives are twice as
#' long as the duration of the available study data, or `2*max(Time_trans)`.
#' Under this assumption, the corresponding elimination and absorption time
#' constants would be `log(2)/(2*max(Time_trans))`. Therefore, the default lower
#' bounds for `kelim` and `kgutabs` are `log(2)/(2*max(Time_trans))`.
#'
#' Upper bounds are based on the opposite assumption: that elimination and
#' absorption are very fast compared to the time scale of the study.
#' Specifically, the upper bounds assume that the elimination and absorption
#' half-lives are half as long as the time of the first observation after time
#' 0, or `0.5*min(Time_trans[Time_trans>0])`. Under this asumption, the
#' corresponding elimination and absorption time constants would be
#' `log(2)/(0.5*min(Time_trans[Time_trans>0]))`. Therefore, the default lower
#' bounds for `kelim` and `kgutabs` are
#' `log(2)/(0.5*min(Time_trans[Time_trans>0]))`.
#' }
#' \subsection{Default lower and upper bounds for `Vdist`}{
#' By default, the lower bound for `Vdist` is 0.01, and the upper bound for
#' `Vdist` is 100. These values were chosen based on professional judgment.
#' }
#' \subsection{Default lower and upper bounds for `Fgutabs`}{
#' By default, the lower bound for `Fgutabs` is 0.0, and the upper bound for
#' `Fgutabs` is 1. These are simply the bounds of the physically-meaningful
#' range for a fraction.
#' }
#' \subsection{Default lower and upper bounds for `Fgutabs_Vdist`}{
#' By default, the lower bound for the ratio `Fgutabs_Vdist` is 0.01, and the
#' upper bound is 100. These values were chosen based on professional judgment.
#' }
#' \subsection{Default lower and upper bounds for `Rblood2plasma`}{
#' By default, the lower bound for the blood:plasma partition coefficient
#' `Rblood2plasma` is 0.01, and the upper bound is 100. These values were chosen
#' based on professional judgment.
#' }
#' \subsection{Default lower and upper bounds for `Frec`}{
#' By default, the lower bound for `Frec` is 0.0001, and the upper bound for `Frec`
#' is 1. This assumes that if `excreta` have any measurable values then there is
#' at the very least there is 0.1% recovery.
#' }
#'
#' @section Starting values for each parameter:
#'
#' Starting values for each parameter (starting guesses for the numerical
#' optimizer) are derived from the data using [get_starts_1comp()].
#'
#' If the starting values returned by [get_starts_1comp()] fall outside the
#' bounds for any parameter(s), then the starting value will be reset to a value
#' halfway between the lower and upper bounds for that parameter.
#'
#' @inheritParams get_params_flat
#' @param restrictive A boolean value (Default: FALSE) that determines whether
#'  to assume restrictive clearance when setting starting values for parameters.
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
#' @author Gilberto Padilla Mercado
#' @family 1-compartment radiation model functions
#' @family get_params functions
#' @family built-in model functions

get_params_1comp_rad <- function(
    data,
    lower_bound = ggplot2::aes(Q_totli = NA,
                               Q_gfr = NA,
                               Q_alv = NA,
                               Kblood2air = NA,
                               Fup = 0,
                               Clint = 0,
                               Vdist = NA,
                               Fgutabs = NA,
                               kgutabs = NA,
                               Rblood2plasma = NA,
                               Frec = 0.001),
    upper_bound = ggplot2::aes(Q_totli = NA,
                               Q_gfr = NA,
                               Q_alv = NA,
                               Kblood2air = NA,
                               Fup = 1,
                               Clint = 1E5,
                               Vdist = NA,
                               Fgutabs = NA,
                               kgutabs = NA,
                               Rblood2plasma = NA,
                               Frec = 1),
    param_units = ggplot2::aes(Q_totli = "L/h/kg ^3/4",
                               Q_gfr = "L/h/kg ^3/4",
                               Q_alv = "L/h/kg ^3/4",
                               Kblood2air = "unitless ratio",
                               Fup = "unitless fraction",
                               Clint = "L/h/kg",
                               Vdist = paste0("(", # Vdist
                                              unique(Dose.Units),
                                              ")",
                                              "/",
                                              "(",
                                              unique(Conc.Units),
                                              ")"),
                               Fgutabs = "unitless fraction", # Fgutabs
                               kgutabs = paste0("1/", # kgutabs
                                                unique(Time_trans.Units)),
                               Rblood2plasma = "unitless ratio",
                               Frec = "unitless fraction"),
    restrictive = TRUE) {
  # param names
  param_name <- c("Q_totli",
                  "Q_gfr",
                  "Q_alv",
                  "Kblood2air",
                  "Fup",
                  "Clint",
                  "Vdist",
                  "Fgutabs",
                  "kgutabs",
                  "Rblood2plasma",
                  "Frec")

  # Default lower bounds, to be used in case the user specified a non-default
  # value for the `lower_bound` argument, but did not specify expressions for all
  # parameters. Any parameters not specified in the `lower_bound` argument will
  # take their default lower bounds defined here. This should be the same as the
  # default value for the `lower_bound` argument.
  lower_bound_default = ggplot2::aes(Q_totli = NA,
                                     Q_gfr = NA,
                                     Fup = 0,
                                     Clint = 0,
                                     Vdist = 0.01,
                                     Fgutabs = 0,
                                     kgutabs = log(2) / (2 * max(Time_trans)),
                                     Fgutabs_Vdist = 0.01,
                                     Rblood2plasma = 1e-2,
                                     Frec = 0.001)
  # which parameters did not have lower bounds specified in the `lower_bound`
  # argument?
  lower_bound_missing <- setdiff(names(lower_bound_default),
                                 names(lower_bound))
  # fill in the default lower bounds for any parameters that don't have them
  # defined in the `lower_bound` argument
  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]

  # Default upper bounds, to be used in case the user specified a non-default
  # value for the `upper_bound` argument, but did not specify expressions for all
  # parameters. Any parameters not specified in the `upper_bound` argument will
  # take their default upper bounds defined here. This should be the same as the
  # default value for the `upper_bound` argument.
  upper_bound_default = ggplot2::aes(Q_totli = NA,
                                     Q_gfr = NA,
                                     Fup = 1,
                                     Clint = 1E5,
                                     Vdist = 100,
                                     Fgutabs = 1,
                                     kgutabs = log(2) / (0.5 * min(Time_trans[Time_trans > 0])),
                                     Fgutabs_Vdist = 1e2,
                                     Rblood2plasma = 100,
                                     Frec = 1)
  # which parameters did not have upper bounds specified in the `upper_bound`
  # argument?
  upper_bound_missing <- setdiff(names(upper_bound_default),
                                 names(upper_bound))
  # fill in the default upper bounds for any parameters that don't have them
  # defined in the `upper_bound` argument
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]

  # initialize optimization: start with optimize = TRUE for all params
  optimize_param <- rep(TRUE, length(param_name))
  # but then will hold constant some of these
  optimize_param[param_name %in% c("Fup", "Q_gfr", "Q_totli", "Rblood2plasma")] <- FALSE

  # initialize whether each param is used: start with use = TRUE for all params
  use_param <- rep(TRUE, length(param_name))

  # now follow the logic described in the documentation for this function:
  if ("oral" %in% data$Route) {
    # if yes oral data:
    if ("iv" %in% data$Route) {
      # if both oral and IV data, Fgutabs and Vdist can be fit separately, so turn off Fgutabs_Vdist
      optimize_param[param_name %in% "Fgutabs_Vdist"] <- FALSE
      use_param[param_name %in% "Fgutabs_Vdist"] <- FALSE
    } else {
      # if oral ONLY:
      # cannot fit Fgutabs and Vdist separately, so turn them off.
      optimize_param[param_name %in% c("Fgutabs", "Vdist")] <- FALSE
      use_param[param_name %in% c("Fgutabs", "Vdist")] <- FALSE
    }
  } else {
    # if no oral data, can't fit kgutabs, Fgutabs, or Fgutabs_Vdist,
    # and they won't be used.
    optimize_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_Vdist")] <- FALSE
    use_param[param_name %in% c("kgutabs", "Fgutabs", "Fgutabs_Vdist")] <- FALSE
  }

  # if both "blood" and "plasma" are not in data,
  # then Rblood2plasma will not be optimized
  if (!(all(c("blood", "plasma") %in% data$Media))) {
    optimize_param[param_name %in% "Rblood2plasma"] <- FALSE
  }

  param_units_vect <- sapply(param_units,
                             function(x) rlang::eval_tidy(x,
                                                          data = data),
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  param_units_vect <- param_units_vect[param_name]

  lower_bound_vect <- sapply(lower_bound,
                             function(x) rlang::eval_tidy(x,
                                                          data = data),
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  lower_bound_vect <- lower_bound_vect[param_name]

  upper_bound_vect <- sapply(upper_bound,
                             function(x) rlang::eval_tidy(x,
                                                          data = data),
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  upper_bound_vect <- upper_bound_vect[param_name]

  par_DF <- data.frame("param_name" = param_name,
                       "param_units" = param_units_vect,
                       "optimize_param" = optimize_param,
                       "use_param" = use_param,
                       "lower_bound" = lower_bound_vect,
                       "upper_bound" = upper_bound_vect)

  # Get starts goes here and the rest of the function
  par_DF <- get_starts_1comp_rad(
    data = data,
    par_DF = par_DF,
    restrictive = restrictive
  )

  # check to ensure starting values are within bounds
  # if not, then replace them by a value halfway between bounds
  # Note here that parameter lower and upper bounds can be EQUAL to start value
  # This is useful for non-optimized parameters
  start_low <- (par_DF$start < par_DF$lower_bound) %in% TRUE
  start_high <- (par_DF$start > par_DF$upper_bound) %in% TRUE
  start_nonfin <- !is.finite(par_DF$start)

  # Notify user which starts are out of bounds:
  if (any(start_low | start_high | start_nonfin)) {

    oob_starts <- par_DF[start_low | start_high,][["param_name"]]
    inf_starts <- par_DF[start_nonfin,][["param_name"]]
    if (!all(start_nonfin)) {
      message("There are out of bounds starting values detected for ",
              paste(unique(data$Chemical), unique(data$Species), collapse = "|"),
              "\n",
              "See summary below:"
      )

      if (length(oob_starts) > 0) {
        message("Out of bounds parameter starts: ",
                paste(oob_starts, collapse = ", "))
      }
      if (length(inf_starts) > 0) {
        message("Non-finite parameter starts: ",
                paste(inf_starts, collapse = ", "))
      }

      message("Values will be replaced by the mean of lower and upper bounds ",
              "when possible.\n"
      )
    } else {
      message(
        "Unable to parameterize ",
        paste(unique(data$Chemical), unique(data$Species), collapse = "|")
      )
    }

  }

  par_DF[start_low | start_high | start_nonfin,
         "start"] <- rowMeans(
           cbind(par_DF[start_low | start_high | start_nonfin,
                        c("lower_bound",
                          "upper_bound")]
           )
         )

  return(par_DF)
}
