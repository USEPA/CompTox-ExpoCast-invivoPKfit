#' Get an elbow point
#'
#' Given a set of data specified as two vectors of `x` and `y` values, find an
#' elbow point.
#'
#' This is a helper function for [get_starts()] to find elbow points.
#'
#' Given a set of (x,y) data points, an "elbow point" can be defined by drawing
#' a line connecting the points for minimum and maximum x, and then finding the
#' x value of the observation where the distance to that line is greatest.
#'
#' [get_starts()] uses elbow points as a way to automate separation of
#' concentration-time data into different kinetic phases in order to calculate
#' starting points for fitting TK model parameters. For example, if
#' concentration-time data are described by a two-compartment TK model, then
#' early and late elimination phases will be separated by an elbow point. This
#' helper function finds the elbow points. (When this function is called from
#' [get_starts()], `x` will be a vector of time values, and `y` will be a vector
#' of log-transformed dose-normalized concentration values.)
#'
#'
#'
#' @param x The `x` values from the data where an elbow point is to be found.
#' @param y The `y` values from the data where an elbow point is to be found.
#' @param ... Optional: additional arguments that will be passed to
#'   [stats::approx()] if it is used.
#' @return A list with two named numeric scalar elements, `x` and `y`. `x`
#'   contains the `x` value at the elbow point. `y` contains the `y` value at
#'   the elbow point. The elbow point is *not* necessarily one of the input data
#'   points; it may be interpolated.
#' @author Caroline Ring
get_elbow <- function(x, y, ...){
  #find elbow point:
  #draw a line connecting first and last observation
  #find the point that has the largest residual
  xmin <- min(x)
  xmax <- max(x)

  obs_first <- c(xmin, median(y[x==xmin]))
  obs_last <- c(xmax, median(y[x==xmax]))
  slope <- (obs_last[2] - obs_first[2])/(obs_last[1] - obs_first[1])
  intercept <- obs_first[2] - slope * obs_first[1]

  med_y <- sapply(x,
                          function(this_x) median(y[x==this_x])
  )

  pred <- slope * x + intercept
  resid <- med_y - pred
  elbow_x <- x[which.max(abs(resid))]

  return(elbow_x)
}

