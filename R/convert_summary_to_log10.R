#' Convert sample mean and SD to log10-scale
#'
#' Estimate log10-scale sample mean and standard deviation from natural-scale
#' sample mean and standard deviation
#'
#' \eqn{\bar{y}_i} is the natural-scale sample mean for group \eqn{i}. \eqn{s_i}
#' is the natural-scale sample standard deviation for group \eqn{i}.
#'
#' \deqn{\textrm{log10-scale sample mean}_i = \log_{10}
#' \left(\frac{\bar{y}_i^2}{\sqrt{\bar{y}_i^2 + s_i^2}} \right)}
#'
#' \deqn{
#' \textrm{log10-scale sample SD}_i =
#'  \sqrt{
#'  \log_{10} \left(1 + \frac{s_i^2}{\bar{y}_i^2} \right)
#'  }
#'  }
#'
#' @param sample_mean Numeric: one or more sample means
#' @param sample_SD Numeric: one or more sample SDs
#' @return A list with two named elements: "log10mean" and "log10SD", the log10-scale
#'   sample means and log10-scale sample SDs, respectively.
#' @export
#' @author Caroline Ring
convert_summary_to_log10 <- function(sample_mean, sample_SD) {
  log10mean <- log10(sample_mean^2 / sqrt(sample_SD^2 + sample_mean^2))
  log10SD <- sqrt(log10(1 + (sample_SD^2 / sample_mean^2)))
  return(list("log10mean" = log10mean, "log10SD" = log10SD))
}
