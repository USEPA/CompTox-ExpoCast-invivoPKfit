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

auto_units <- function(y,
                       from){
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

    target <- 10
    #auto-select units based on the midpoint of x
    #goal: to keep the midpoint of x somewhere around 10-ish
    #on the grounds that the midpoint of a PK study might be around 1 elimination halflife
    #and absorption halflife might be around 1/5 of elimination halflife
    #so, midpoint of 10 puts both half-lives on a potentially reasonable scale

    y_mid <- mean(range(y))
    #try conversion to each of the period_units
    #starting with the ones above and below the "from" units --
    #go in the direction of the one that gets closer to putting the midpoint at 5 --
    #then continue in that direction until we start getting farther away from 5 again

    i <- which(period_units %in% from)
    target_dist <- abs(y_mid - target)


    #first we choose a direction -- up or down
    if(i == 1){
      #if we're at the low end we can only go up
      dir <- 1
      y_mid_up <- mean(
        range(
          as.numeric(
            lubridate::duration(num = y,
                                units = from),
            period_units[i + 1])
        )
      )
      target_dist_new <- abs(y_mid_up - target)
    }else if(i == length(period_units)){
      #if we're at the high end we can only go down
      dir <- -1
      y_mid_down <- mean(
        range(
          as.numeric(
            lubridate::duration(num = y,
                                units = from),
            period_units[i - 1])
        )
      )

      target_dist_new <- abs(y_mid_down - target)
    }else{
      #try both ways and see which one gets closer to a midpoint of target
      y_mid_up <- mean(
        range(
          as.numeric(
            lubridate::duration(num = y,
                                units = from),
            period_units[i + 1])
        )
      )
      target_dist_up <- abs(y_mid_up - target)

      y_mid_down <- mean(
        range(
          as.numeric(
            lubridate::duration(num = y,
                                units = from),
            period_units[i - 1])
        )
      )

      target_dist_down <- abs(y_mid_down - target)

      if(target_dist_up < target_dist_down){
        #if up gets closer to target than down does
        dir <- 1
        target_dist_new <- target_dist_up
      }else if(target_dist_down < target_dist_up){
        #if down gets closer to target than up does
        dir <- -1
        target_dist_new <- target_dist_down
      }
    }

    #now search in that direction until either we hit the end, or until it starts getting *farther* away from target
    i <- i + dir
    #are we getting closer (negative delta) or farther away (positive delta)?
    delta_target_dist <- target_dist_new - target_dist
    target_dist <- target_dist_new

    while(i > 1 &
          i < length(period_units) &
          delta_target_dist < 0){
      #try the next step
      i_new <- i + dir
      y_mid_new <- mean(
        range(
          as.numeric(
            lubridate::duration(num = y,
                                units = from),
            period_units[i_new])
        )
      )
      target_dist_new <- abs(y_mid_new - target)
      #are we getting closer (negative delta) or farther away (positive delta)?
      delta_target_dist <- target_dist_new - target_dist

      i <- i_new
      target_dist <- target_dist_new

    }

    #select the best unit
    to <- period_units[i - dir]

  return(to)
}

