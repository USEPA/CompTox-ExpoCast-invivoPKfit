#' Scale concentrations
#'
#'
#' @param dose_norm Logical: Whether to normalize observed concentrations (and
#'  observed concentration standard deviations and limits of quantification) by
#'  dividing them by the corresponding dose. Default `FALSE`.
#' @param log10_trans Logical: Whether to apply a `log10()` transformation to
#'  observed concentrations (and limits of quantification), after any dose
#'  normalization is applied. Default `FALSE`.
#' @param ... Other arguments (not currently used)
#' @return An object of class `pk_scales`: A named list with elements supplied to
#'  [scale_conc()]). This object is usually added to an existing [pk()] object
#'  using `+`. See [pk_add.pk_scales()].
#' @export
#' @author Caroline Ring
scale_conc <- function(dose_norm = FALSE, log10_trans = FALSE, ...) {

  expr <- quote(.conc)
  if (dose_norm %in% TRUE) expr <- substitute(x / Dose, env = list(x = expr))
  if (log10_trans %in% TRUE) expr <- substitute(log10(x), env = list(x = expr))

  # Initialize scale_conc as a list with element "name"
  scale_conc <- list(name = "conc",
                     value = list(
                     "dose_norm" = dose_norm,
                     "log10_trans" = log10_trans,
                     "expr" = expr)
  )

  # get any other arguments and values
  scale_conc$value <- c(scale_conc$value, list(...))
  # set class
  class(scale_conc) <- c(class(scale_conc), "pkproto", "pk_scales")

return(scale_conc)
}

#' Scale times
#'
#' Transform time data
#'
#' @param new_units New units to use for time. Default is `"identity"` (leave
#'   time in the original units). Another useful option is `"auto"`, to
#'   automatically select new time units based on the time of the last detected
#'   observation. You may also specify any time units understood by
#'   `lubridate::duration()`, i.e., `"seconds"`, `"hours"`, `"days"`, `"weeks"`,
#'   `"months"`, `"years"`, `"milliseconds"`, `"microseconds"`, `"nanoseconds"`,
#'   and/or `"picoseconds"`. You may only specify one new unit (e.g., `new_units
#'   = c("days", "weeks")` is not valid).
#' @param ... Other arguments (not currently used)
#' @return An object of class `pk_scales`: A named list with two elements `name
#'   = "time"` (denoting the variable to be scaled) and `value = list("new_units" =
#'   new_units, ...)` (denoting the arguments supplied to [scale_time()]). See
#'   [pk_add.pk_scales()].
scale_time <- function(new_units = "identity",
                          ...) {
  # get arguments and values
  argg <- c(as.list(environment()), list(...))
  scale_time <- list(name = "time")
  scale_time$value <- argg
  # set class
  class(scale_time) <- c(class(scale_time), "pkproto", "pk_scales")

  return(scale_time)
}
