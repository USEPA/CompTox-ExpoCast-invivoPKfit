#' Helper function to find peak
#'
#' Finds x- and y-value at peak y value.
#'
#' This is a helper function called by [get_starts()] to estimate the end of the
#' absorption phase for oral TK data.
#'
#' If there is more than one unique `x` value where both `x` and corresponding
#' `y` are finite, this function calls [stats::approx()], then uses
#' [base::which.max()] to locate the maximum interpolated `y`-value. In effect,
#' this takes the average `y` value for each unique `x` value, then finds the
#' maximum.
#'
#' If there is only one unique `x` value where both `x` and corresponding `y`
#' are finite, this function returns that `x` value, along with the average of
#' the finite `y` values corresponding to it.
#'
#' If there are no unique `x` values where both `x` and corresponding `y` are
#' finite, this function returns `NA_real_` for the peak `x` and `y` values.
#'
#' @param x A numeric vector of `x` data
#' @param y A numeric vector of `y` data
#' @param ... Optional: Additional arguments which will be passed to
#'   [stats::approx()] (other than `x`, `y`, and `xout`).
#' @return A list with two named numeric scalar components, `x` and `y`,
#'   containing the x- and y-values at the peak.
#' @author Caroline Ring

get_peak <- function(x, y, ...){

  x_finite <- x[is.finite(x) & is.finite(y)]
  y_finite <- y[is.finite(x) & is.finite(y)]

  if(length(unique(x_finite))>1){
  suppressWarnings(tmp <- approx(x=x, y=y, xout = unique(x), ...))
    peak_ind <- which.max(tmp$y)

    peak <- list("x" = tmp$x[peak_ind],
                "y" = tmp$y[peak_ind])

  }else if(length(unique(x_finite))==1){
    peak <- list("x" = unique(x_finite),
                "y" = mean(y_finite))
  }else{
    peak <- list("x" = NA_real_,
                "y" = NA_real_)
  }

return(peak)
}
