#'Automatically select new time units
#'
#'Given a vector of time values in original units, this function selects new
#'time units such that, when time is rescaled to the new units, the midpoint of
#'the time vector is as close to 10 as possible.
#'
#' # Acceptable/understood time units in `period_units`
#'
#' ```
#' c("picoseconds",
#'"nanoseconds",
#'"microseconds",
#'"milliseconds",
#'"seconds",
#'"minutes",
#'"hours",
#'"days",
#'"weeks",
#'"months",
#'"years")
#' ```
#'
#'@param y A numeric vector of time values
#'@param from The original units of `y`
#'@param target The target value (order of magnitude) for the midpoint of rescaled time values. Default 10.
#'@param period_units A list of acceptable/understood time units. See Details. Default `time_units`.
#'@return Character: Automatically-selected new time units, which will be one of `period_units`.
#'@export
#'@author Caroline Ring
auto_units <- function(y,
                       from,
                       target = 10,
                       period_units = time_units){


  #auto-select units based on the midpoint of x
  #goal: to keep the midpoint of x somewhere around 10-ish
  #on the grounds that the midpoint of a PK study might be around 1 elimination halflife
  #and absorption halflife might be around 1/5 of elimination halflife
  #so, midpoint of 10 puts both half-lives on a potentially reasonable scale

  y_mid <- midpt_log10(y)
  target_log10 <- log10(target)
  #try conversion to each of the period_units
  #starting with the ones above and below the "from" units --
  #go in the direction of the one that gets closer to putting the midpoint at 5 --
  #then continue in that direction until we start getting farther away from 5 again

  i <- which(period_units %in% from)

  target_dist <- abs(y_mid - target_log10)

  #first we choose a direction -- up or down
  if (i == 1) {
    #if we're at the low end we can only go up
    dir <- 1
    y_mid_up <- midpt_log10(
      convert_time(y,
                   from = from,
                   to = period_units[i+1])
    )

    target_dist_new <- abs(y_mid_up - target_log10)
  } else if (i == length(period_units)) {
    #if we're at the high end we can only go down
    dir <- -1
    y_mid_down <- midpt_log10(
      convert_time(y,
                   from = from,
                   to = period_units[i-1])
    )

    target_dist_new <- abs(y_mid_down - target_log10)
  } else {
    #try both ways and see which one gets closer to a midpoint of target
    y_mid_up <- midpt_log10(
      convert_time(y,
                   from = from,
                   to = period_units[i+1])
    )
    target_dist_up <- abs(y_mid_up - target_log10)

    y_mid_down <- midpt_log10(
      convert_time(y,
                   from = from,
                   to = period_units[i-1])
    )

    target_dist_down <- abs(y_mid_down - target_log10)

    if (target_dist_up < target_dist_down) {
      #if up gets closer to target than down does
      dir <- 1
      target_dist_new <- target_dist_up
    } else if (target_dist_down < target_dist_up) {
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

  while (i > 1 & i < length(period_units) & delta_target_dist < 0) {
    #try the next step
    i_new <- i + dir
    y_mid_new <- midpt_log10(
      convert_time(y,
                   from = from,
                   to = period_units[i_new])
    )
    target_dist_new <- abs(y_mid_new - target_log10)
    #are we getting closer (negative delta) or farther away (positive delta)?
    delta_target_dist <- target_dist_new - target_dist

    i <- i_new
    target_dist <- target_dist_new

  }

  #select the best unit
  to <- period_units[i - dir]

  return(to)
}

#' Log10-scaled midpoint
#'
#' @param x A numeric vector
#' @return The log10-scaled midpoint, calculated as `log10(mean(range(x))`
midpt_log10 <- function(x){
  log10(mean(range(x)))
}

