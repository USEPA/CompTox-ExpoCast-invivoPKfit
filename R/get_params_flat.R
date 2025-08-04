#' Get parameters to be optimized for flat model
#'
#' The full set of model parameters for the flat model includes
#' `Vdist`, `Fgutabs`, and `Rblood2plasma`. Whether each one can be estimated
#' from the data depends on what routes of administration are included in the
#' data.
#'
#' @section IV data, no oral data:
#'
#' If IV dosing data are available, but no oral dosing data are available, then
#' only the parameter `Vdist` will be estimated from the data. The
#' parameter `Fgutabs` cannot be estimated from IV data alone and will not be used to evaluate the model.
#'
#' @section Oral data, no IV data:
#'
#' If oral dosing data are available, but no IV dosing data are available, then
#' the parameters `Fgutabs` and `Vdist` cannot be identified separately. From
#' oral data alone, only the ratio `Fgutabs/Vdist` can be identified. This ratio
#' is represented by a single parameter named `Fgutabs_Vdist`. `Fgutabs` and
#' `Vdist` will not be estimated nor used in model evaluation, but `Fgutabs_Vdist` will be estimated.
#'
#' @section Oral data and IV data:
#'
#' If both oral and IV dosing data are available, then `Vdist` and `Fgutabs` will both be estimated from the data.
#'
#' @section Blood and plasma data:
#'
#' If both blood and plasma data are available, then `Rblood2plasma` will be estimated from the data.
#'
#' @section Only one of blood or plasma data:
#'
#' If only one of blood or plasma data are available, then `Rblood2plasma` will be
#' held constant at 1, not estimated from the data.
#'
#' @section Default lower and upper bounds for each parameter:
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
#'
#'
#' @section Starting values for each parameter:
#'
#' Starting values for each parameter (starting guesses for the numerical
#' optimizer) are derived from the data using [get_starts_flat()].
#'
#' If the starting values returned by [get_starts_flat()] fall outside the
#' bounds for any parameter(s), then the starting value will be reset to a value
#' halfway between the lower and upper bounds for that parameter.
#'
#' @param data The data set to be fitted (e.g. the result of [preprocess_data()])
#' @param lower_bound A mapping specified using a call to [alist()],
#'  giving the lower bounds for each variable, as expressions which may include
#'  variables in `data`.
#' @param upper_bound A mapping specified using a call to [alist()],
#'  giving the upper bounds for each variable, as expressions which may include
#'  variables in `data`.
#' @param param_units A mapping specified using a call to [alist()],
#'  giving the units for each variable, as expressions which may include
#'  variables in `data`.
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
#' @family flat model functions
#' @family get_params functions
#' @family built-in model functions

get_params_flat <- function(
    data,
    lower_bound = NULL,
    upper_bound = NULL,
    param_units = alist(
      Vdist = paste0("(", unique(Dose.Units), ")/(", unique(Conc.Units), ")"),
      Fgutabs = "unitless fraction",
      Fgutabs_Vdist = paste0("(", unique(Conc.Units), ")/(", unique(Dose.Units), ")"),
      Rblood2plasma = "unitless ratio"
    )
) {
  # param names
  param_name <- c("Vdist", "Fgutabs", "Fgutabs_Vdist", "Rblood2plasma")

  lower_bound_default = alist(
    Vdist = 0.01,
    Fgutabs = 1E-3,
    Fgutabs_Vdist = 1E-5,
    Rblood2plasma = 1e-2
  )

  lower_bound_missing <- base::setdiff(names(lower_bound_default),
                                 names(lower_bound))
  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]

  upper_bound_default = alist(
    Vdist = 100,
    Fgutabs = 1,
    Fgutabs_Vdist = 1e2,
    Rblood2plasma = 100
  )

  upper_bound_missing <- base::setdiff(names(upper_bound_default),
                                 names(upper_bound))
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]

  # Set the parameters to optimize and use here
  optimize_param <- rep(TRUE, length(param_name))
  use_param <- rep(TRUE, length(param_name))

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
  # if no oral data, can't fit Fgutabs, or Fgutabs_Vdist,
  # and they won't be used.
  optimize_param[param_name %in% c("Fgutabs", "Fgutabs_Vdist")] <- FALSE
  use_param[param_name %in% c("Fgutabs", "Fgutabs_Vdist")] <- FALSE
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
                             USE.NAMES = TRUE)
  param_units_vect <- param_units_vect[param_name]

  lower_bound_vect <- sapply(lower_bound,
                             rlang::eval_tidy,
                             data = data,
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  lower_bound_vect <- lower_bound_vect[param_name]

  upper_bound_vect <- sapply(upper_bound,
                             rlang::eval_tidy,
                             data = data,
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  upper_bound_vect <- upper_bound_vect[param_name]

  par_DF <- data.frame("param_name" = param_name,
                       "param_units" = param_units_vect,
                       "optimize_param" = optimize_param,
                       "use_param" = use_param,
                       "lower_bound" = lower_bound_vect,
                       "upper_bound" = upper_bound_vect)

  rownames(par_DF) <- par_DF$param_name

  par_DF <- get_starts_flat(data = data,
                            par_DF = par_DF)

  # check to ensure starting values are within bounds
  # if not, then replace them by a value halfway between bounds
  start_low <- (par_DF$start < par_DF$lower_bound) %in% TRUE
  start_high <- (par_DF$start > par_DF$upper_bound) %in% TRUE
  start_nonfin <- !is.finite(par_DF$start)

  par_DF[start_low | start_high | start_nonfin, "start"] <- rowMeans(
    cbind(par_DF[start_low | start_high | start_nonfin,
                 c("lower_bound",
                   "upper_bound")]
    )
  )

  par_DF$start <- as.list(par_DF$start)

  return(par_DF)

}
