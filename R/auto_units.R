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

