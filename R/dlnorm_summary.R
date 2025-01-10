#' Log-normal distribution density function for summary data
#'
#' Evaluates the normal distribution density function for summary data reported
#' as sample mean, sample SD, and sample N. Sample mean and sample SD should be
#' on the *natural* scale. If you have log-scale sample mean and SD (i.e., the
#' mean and SD of log-transformed observations),then use [dnorm_summary()]
#' instead.
#'
#' `x_mean`, `x_sd`, `X_N`, `mu`, and `sigma` should either be all the same
#' size, or length 1. If they are different lengths, they will be repeated until
#' their lengths match, with a warning.
#'
#'
#' @param mu *Log-scale* mean of the log-normal distribution to be evaluated
#'   (*not* the sample mean). May be a numeric scalar or vector.
#' @param sigma *Log-scale* standard deviation of the log-normal distribution to
#'   be evaluated (*not* the sample SD). May be a numeric scalar or vector.
#' @param x_mean Sample mean (on the *natural* scale). May be a numeric scalar
#'   or vector.
#' @param x_sd Sample standard deviation (on the *natural* scale). May be a
#'   numeric scalar or vector.
#' @param x_N Sample number of observations. May be a numeric scalar or vector.
#' @param log TRUE/FALSE: Whether to return the log of the density function
#'   (i.e., the log-likelihood). Default FALSE.
#'
#' @return A numeric scalar or vector matching the length of the longest of
#'   `mu`, `sigma`, `x_mean`, `x_sd`, and `x_N`.
#' @author Caroline Ring
#' @export

dlnorm_summary <- function(mu,
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
            paste(paste0(names(x_len)[bad_len], " length = ", x_len[bad_len]),
                  collapse = "\n"
            ),
            "\n They will be repeated to match the length of the longest input,",
            names(x_len)[which_max_len],
            " length = ",
            max_len,
            "."
    )
    # repeat to match longest
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
    + (-1 / (2 * sigma^2))
    * ((x_N - 1) * log(1 + x_sd^2 / x_mean^2)
       + x_N * (log(x_mean^2 / (sqrt(x_sd^2 + x_mean^2))))^2)
    + mu / sigma^2 * x_N * log(x_mean^2 / (sqrt(x_sd^2 + x_mean^2)))
    - mu^2 * x_N / (2 * sigma^2)
  )

 if (log == TRUE) {
   return(y_log)
 } else {
   return(exp(y_log))
 }

}
