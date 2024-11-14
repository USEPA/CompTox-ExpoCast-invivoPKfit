#' Find the peak of a data series
#'
#' Finds x- and y-value at peak y value.
#'
#' If there is more than one unique `x` value where both `x` and corresponding
#' `y` are finite, this function calls [stats::approx()] with `method = 'linear'`, then uses
#' [base::which.max()] to locate the maximum interpolated `y`-value.
#'
#' If there is only one unique `x` value where both `x` and corresponding `y`
#' are finite, this function calls [stats::approx()] with `method = 'constant'`, then uses
#' [base::which.max()] to locate the maximum interpolated `y`-value.
#'
#' If there are no unique `x` values where both `x` and corresponding `y` are
#' finite, this function returns `NA_real_` for the peak `x` and `y` values.
#'
#' @param x A numeric vector of `x` data
#' @param y A numeric vector of `y` data
#' @param ties As for [stats::approxfun()]: The function to apply to y-values
#'   that have the same x-value. Default `'median'`. `'mean'` may also be
#'   useful.
#' @param na.rm As for [stats::approxfun()]: How to handle missing values.
#'   Default `TRUE` to exclude missing values from analysis.
#' @param ... Optional: Additional arguments which will be passed to
#'   [stats::approx()] (other than `x`, `y`, and `xout`).
#' @return A list with two named numeric scalar components, `x` and `y`,
#'   containing the x- and y-values at the peak.
#' @author Caroline Ring

get_peak <- function(x, y, ties = "median", na.rm = TRUE, ...) {

  x_finite <- x[is.finite(x) & is.finite(y)]
  y_finite <- y[is.finite(x) & is.finite(y)] # Not used?

  if (length(unique(x_finite)) > 1) {
    tmp <- suppressWarnings(approx(x = x,
                                   y = y,
                                   xout = unique(x),
                                   method = "linear",
                                   ties = ties,
                                   na.rm = na.rm,
                                   ...))
    peak_ind <- which.max(tmp$y)

    peak <- list("x" = tmp$x[peak_ind], "y" = tmp$y[peak_ind])

  } else if (length(unique(x_finite)) == 1) {
    tmp <- suppressWarnings(approx(x = x,
                                   y = y,
                                   xout = unique(x),
                                   method = "constant",
                                   ties = ties,
                                   na.rm = na.rm,
                                   ...))
    peak_ind <- which.max(tmp$y)

    peak <- list("x" = tmp$x[peak_ind],
                 "y" = tmp$y[peak_ind])

  } else {
    peak <- list("x" = NA_real_,
                 "y" = NA_real_)
  }

  return(peak)
}
