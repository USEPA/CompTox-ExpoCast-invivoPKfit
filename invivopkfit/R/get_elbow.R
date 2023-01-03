#' Get an elbow point
#'
#' Given a set of data specified as two vectors of `x` and `y` values, find an
#' elbow point.
#'
#' This is a helper function for `get_starts()` to find elbow points.
#'
#' Given a set of (x,y) data points, an "elbow point" is qualitatively defined
#' as the (x,y) point where the slope changes.
#'
#' `get_starts()` uses elbow points as a way to automate separation of
#' concentration-time data into different kinetic phases in order to calculate
#' starting points for fitting TK model parameters. For example, if
#' concentration-time data are described by a two-compartment TK model, then
#' early and late elimination phases will be separated by an elbow point. This
#' helper function finds the elbow points. (When this function is called from
#' `get_starts()`, `x` will be a vector of time values, and `y` will be a vector
#' of log-transformed dose-normalized concentration values.)
#'
#' If there are at least four unique `x` values whose corresponding `y` values
#' are finite, `get_elbow()` calls `akmedoids::elbow_point()`.
#'
#' If there are fewer than four but more than one unique `x` values whose
#' corresponding `y` values are finite, or if `akmedoids::elbow_point()` fails
#' with an error or returns an elbow point whose `x` or `y` value is non-finite,
#' then the elbow point `x` value will be selected using the helper function
#' `get_middle()` as the middle `x` point (halfway between the minimum and
#' maximum `x` values), and the elbow `y`-value will be interpolated at that
#' middle `x` point using `stats::approx()`. This is an arbitrary choice, not
#' meaningful if the goal is to locate an actual elbow point, but a useful
#' default when the goal is to split TK data into two phases.
#'
#' If there is only one unique `x` value whose corresponding `y` values are
#' finite, then the elbow point `x` value will be selected as that one unique
#' `x` value, and the elbow point `y` value will be selected as the median of
#' the finite `y` values corresponding to that unique `x` value. Again, this is
#' an arbitrary choice, not meaningful if the goal is to locate an actual elbow
#' point, but a useful default when the goal is to split TK data into two
#' phases.
#'
#' @param x The `x` values from the data where an elbow point is to be found.
#' @param y The `y` values from the data where an elbow point is to be found.
#' @param ... Optional: additional arguments that will be passed to
#'   `stats::approx()` if it is used.
#' @return A list with two components, `x` and `y`. `x` contains the `x` value
#'   at the elbow point. `y` contains the `y` value at the elbow point. The
#'   elbow point is *not* necessarily one of the input data points; it may be
#'   interpolated.

get_elbow <- function(x, y, ...){
  #get points where both x and y are finite
  good_ind <- is.finite(x) & is.finite(y)
  x_good <- x[good_ind]
  y_good <- y[good_ind]

  if(length(unique(x_good))>=4){
    elbow <- tryCatch(akmedoids::elbow_point(x_good,
                                             y_good)[c("x", "y")],
                      error = function(x, y, ...){
                        return(get_middle(x, y, ...))
                      })

    #in case we get NA, NaN, or infinite elbow points,
    #fallback to middle x
    if(any(sapply(elbow, function(z)!is.finite(z)))){
    elbow <- get_middle(x,y, ...)
    }
  }else if(length(unique(x_good))>=2){
    #if fewer than 4, but more than 1 unique x values with finite y values
    elbow <- get_middle(x,y, ...)
  }else if(length(unique(x_good))==1){
    #if only 1 unique x value with finite y value
    elbow <- list(x = unique(x_good),
                  y = median(y_good))
  }else{
#if no x-values with finite y-values, return NA for elbow point
    elbow <- list(x = NA_real_,
                  y = NA_real_)
  }

  return(elbow)
}

#' Helper function to get middle point
#'
#' Given a set of data provided as `x` and `y` vectors, return the middle `x`
#' value and an approximation of the corresponding `y` value at that `x` value.
#'
#' The middle `x` value is defined as `mean(range(x, na.rm = TRUE))`. To
#' approximate the corresponding `y` value, `stats::approx()` is used.
#'
#' This helper function is called by `get_elbow()` in case
#' `akmedoids::elbow_point()` fails.
#'
#' @param x A numeric vector of `x` values
#' @param y A numeric vector of `y` values
#' @param ... Optional additional arguments to `stats::approx()` (other than
#'   `x`, `y`, and `xout`). See documentation for `stats::approx()` for details
#'   and options.
#' @return A list with two numeric scalar components, `x` and `y`. `x` contains
#'   the middle `x` value from the input data. `y` contains the interpolated `y`
#'   value at that middle `x` value.
get_middle <- function(x,
                       y,
                       ...){
  x[!is.finite(x)] <- NA_real_
  y[!is.finite(y)] <- NA_real_
  elbow_x <- mean(range(x, na.rm = TRUE), na.rm = TRUE)
  suppressWarnings(elbow_y <- approx(x = x,
                    y = y,
                    xout = elbow_x,
                    ...)$y)
  return(list(x = elbow_x,
              y = elbow_y))
}
