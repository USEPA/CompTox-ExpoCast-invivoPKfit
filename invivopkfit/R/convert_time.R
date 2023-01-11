#' Helper function to convert time units
#'
#' Convert a vector of times between units
#'
#' This helper function uses `lubridate::duration()` to convert between time
#' units.
#'
#' @param x Numeric: one or more time values to be converted.
#' @param from Character vector: `x` is currently in these units. Must be units
#'   understood by `lubridate::duration()`, i.e. "seconds", "hours", "days",
#'   "weeks", "months", "years", "milliseconds", "microseconds", "nanoseconds",
#'   and/or "picoseconds". Default value is "hours".
#' @param to Character vector: `x` will be converted to these units. Must be
#'   units understood by `lubridate::duration()`, i.e. "seconds", "hours",
#'   "days", "weeks", "months", "years", "milliseconds", "microseconds",
#'   "nanoseconds", and/or "picoseconds". Default value is "hours".
#' @param inverse Logical: TRUE if `x` is in units of *inverse* time (e.g.
#'   1/hour, 1/day); FALSE if `x` is in units of time (e.g. hours, days).
#'   Default value is FALSE.
#' @return A numeric vector the same length as `x`, converted from the units in
#'   `from` to the units in `to`.
#'

convert_time <- function(x,
                                   from = "hours",
                                   to = "hours",
                                   inverse = FALSE){

  if(inverse %in% TRUE){
    y <- 1/x
  }else{y <- x}

  y_new <- as.numeric(
   lubridate::duration(num = y,
                       units = from),
   to)

  if(inverse %in% TRUE){
    yout <- 1/y_new
  }else{
    yout <- y_new
  }

  return(yout)
}

