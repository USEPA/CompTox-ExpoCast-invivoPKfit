match_len <- function(time, dose, route, medium){
  #check whether lengths of time, dose, and route match
  time_len <- length(time)
  dose_len <- length(dose)
  route_len <- length(route)
  medium_len <- length(medium)

  len_all <- c(time_len, dose_len, route_len, medium_len)
  #Cases:
  # All three lengths are the same -- OK
  # Two lengths are the same and the third is 1 -- OK
  # Two lengths are 1 and the third is not 1 -- OK
  # Otherwise there is a problem
  good_len <- (length(unique(len_all)) == 1) |
    (length(unique(len_all)) == 2 &
       sum(len_all == 1) %in% c(1, 2))

  if(!good_len){
    stop(paste0("invivopkfit::match_len(): ",
                "'time', 'dose', 'route', and 'medium'",
                "must either be the same length or length 1.\n",
                "'time' is length ", time_len, "\n",
                "'dose' is length ", dose_len, "\n",
                "'route' is length ", route_len, "\n",
                "'medium' is length ", medium_len, "\n")
    )
  }

  #if any are length-1, repeat them to match the longest
  max_len <- max(len_all)
  time <- rep(time, length.out = max_len)
  dose <- rep(dose, length.out = max_len)
  route <- rep(route, length.out = max_len)
  medium <- rep(medium, length.out = max_len)

  return(list("time" = time,
              "dose" = dose,
              "route" = route,
              "medium" = medium))

}
