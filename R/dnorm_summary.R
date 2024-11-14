#' Normal distribution density function for summary data
#'
#' Evaluates the normal distribution density function for summary data reported
#' as sample mean, sample SD, and sample N.
#'
#' `x_mean`, `x_sd`, `X_N`, `mu`, and `sigma` should either be all the same size, or length 1.
#' If they are different lengths, they will be repeated until their lengths
#' match, with a warning.
#'
#'
#' @param mu Mean of the normal distribution to be evaluated (*not* the sample
#'   mean). May be a numeric scalar or vector.
#' @param sigma Standard deviation of the normal distribution to be evaluated
#'   (*not* the sample SD). May be a numeric scalar or vector.
#' @param x_mean Sample mean. May be a numeric scalar or vector.
#' @param x_sd Sample standard deviation. May be a numeric scalar or vector.
#' @param x_N Sample number of observations. May be a numeric scalar or vector.
#' @param log TRUE/FALSE: Whether to return the log of the density function.
#'   Default FALSE (to return the density function value on the natural scale).
#'
#' @return A numeric scalar or vector matching the length of the longest of
#'   `mu`, `sigma`, `x_mean`, `x_sd`, and `x_N`.
#' @author Caroline Ring
#' @export

dnorm_summary <- function(mu,
                             sigma,
                             x_mean,
                             x_sd,
                             x_N,
                             log = FALSE) {

  x_len <- c("x_mean" = length(x_mean),
             "x_sd" = length(x_sd),
             "x_N" = length(x_N),
             "mu" = length(mu),
             "sigma" = length(sigma))

  if (any(x_len %in% 0)) {
    stop("invivopkfit::dnorm_summary(): ",
         "the following arguments have zero length: ",
         toString(names(x_len)[x_len %in% 0])
    )
  }

  max_len <- max(x_len)
  which_max_len <- which.max(x_len)

  bad_len <- (x_len < max_len) & (x_len != 1)


  if (any(bad_len)) {
    warning("invivopkfit::dnorm_summary(): ",
            "the following inputs do not have matching lengths: ",
            paste(paste0(names(x_len)[bad_len],
                         " length = ",
                         x_len[bad_len]),
                  collapse = "\n"
            ),
            "\n They will be repeated to match the length of the longest input,",
            names(x_len)[which_max_len],
            " length = ",
            max_len,
            "."
    )

    # repeat to match longest IFF there is mismatch
    for (i in seq_along(x_len)) {
      assign(names(x_len)[i],
             rep( # repeat the current value of each item to match the length
               get(names(x_len)[i]), # get the current value of each item
               length.out = max_len)
      )
    }
  }

  # Evaluate
  y_log <- (
    x_N * log(1 / (sigma * sqrt(2 * pi)))
    + (-1 / (2 * sigma^2)) * ((x_N - 1) * x_sd^2 + x_N * (x_mean - mu)^2)
  )

 if (log == TRUE) {
   return(y_log)
 } else {
   return(exp(y_log))
 }

}
