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
#'   either "auto", or units understood by `lubridate::duration()`, i.e.
#'   "seconds", "hours", "days", "weeks", "months", "years", "milliseconds",
#'   "microseconds", "nanoseconds", and/or "picoseconds". Default value is
#'   "hours". If "auto", then units will be automatically chosen that make the
#'   midpoint of `x` (or its inverse, if `inverse = TRUE`) as close to 10 as
#'   possible.
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

  period_units <- c("picoseconds",
                    "nanoseconds",
                    "microseconds",
                    "milliseconds",
                    "seconds",
                    "minutes",
                    "hours",
                    "days",
                    "weeks",
                    "months",
                    "years")

  if(inverse %in% TRUE){
    y <- 1/x
  }else{y <- x}


  if(to %in% "auto"){
    #select units automatically
    to <- auto_units(y = y,
                     from = from)
  }


  if(to %in% period_units){
  y_new <- as.numeric(
   lubridate::duration(num = y,
                       units = from),
   to)
  }else{
    if(to %in% "decades"){
      y_new <- as.numeric(
        lubridate::duration(num = y,
                            units = from),
        "years")/10
    }else if(to %in% "centuries"){
      y_new <- as.numeric(
        lubridate::duration(num = y,
                            units = from),
        "years")/100
    }else{
      y_new <- NA_real_
      message(paste("invivopkfit::convert_time(): Allowable values for 'to' are",
                    paste(c(period_units,
                            "decades",
                            "centuries"),
                          collapse = ", "),
                    ", but 'to' =", to,
                    ". Returning NA_real_ for converted times."))
    }
  }

  if(inverse %in% TRUE){
    yout <- 1/y_new
  }else{
    yout <- y_new
  }

  return(yout)
}

