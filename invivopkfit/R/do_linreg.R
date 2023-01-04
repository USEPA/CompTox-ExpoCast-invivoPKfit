#' Fit a line to data
#'
#' Helper function to do linear regression for method of residuals
#'
#' If there are more than two `x` values, this function calls `stats::lm()` and
#' takes the coefficients. If there are exactly two `x` values, this function d
#'
#' @param x A numeric vector of `x` values
#' @param y A numeric vector of `y` values
#' @param intercept_exp Logical: Whether to exp-transform the intercept before
#'   returning. Default TRUE. For example, this is useful when `y` is
#'   log-transformed and you ultimately want the intercept on the natural scale.
#' @param slope_neg Logical: Whether to take the negative of the slope before
#'   returning. Default TRUE. For example, this is useful when a time constant
#'   is the negative of the slope, and you ultimately desire the time constant.
#' @return A list with two numeric scalar components: `intercept` and `slope`,
#'   containing the intercept and slope of the fitted line, respectively
#' @author Caroline Ring
#'
do_linreg <- function(x, y,
                      intercept_exp = TRUE,
                      slope_neg = TRUE){
  x_finite <- x[is.finite(x) & is.finite(y)]
  y_finite <- y[is.finite(x) & is.finite(y)]
  if(length(unique(x_finite))>=2){
    #if there is enough data to do linear regression, do so
    lm_out <- tryCatch(lm(y ~ x),
                       error = function(err) NA_real_)
    if("lm" %in% class(lm_out)){
      intercept <- coef(lm_out)[1]
      slope <- coef(lm_out)[2]
    }else{
      intercept <- NA_real_
      slope <- NA_real_
    }
  }else{
    #if there is only 1 unique x value, or none
    intercept <- NA_real_
    slope <- NA_real_
  }

  if(intercept_exp %in% TRUE){
    intercept <- exp(intercept)
  }

  if(slope_neg %in% TRUE){
    slope <- -1 * slope
  }

  return(list("intercept" = intercept,
              "slope" = slope))
}
