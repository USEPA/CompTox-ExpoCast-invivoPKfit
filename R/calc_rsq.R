#' Calculate r-squared for observed vs. predicted values
#'
#' Calculate the square of the Pearson correlation coefficient (r) between
#' observed and model-predicted values
#'
#' Calculate the square of the Pearson correlation coefficient (r) between
#' observed and model-predicted values, when observed data may be left-censored
#' (non-detect) or may be reported in summary form (as sample mean, sample
#' standard deviation, and sample number of subjects). Additionally, handle the
#' situation when observed data and predictions need to be log-transformed
#' before RMSE is calculated.
#'
#' \eqn{r^2} is calculated according to the following formula, to properly
#' handle observations reported in summary format:
#'
#' \deqn{ r^2 = \left(
#' \frac{
#' \sum_{i=1}^G \mu_i n_i \bar{y}_i -
#' (\bar{\mu} + \bar{y}) \sum_{i=1}^G n_i \mu_i +
#' (\bar{mu} \bar{y}) \sum_{i=1}^G n_i
#' }
#' { \sqrt{
#' \sum_{i=1}^G (n_i - 1) s_i^2 +
#' \sum_{i=1}^G n_i \bar{y}_i^2 -
#' 2 \bar{y} \sum_{i=1}^G n_i \bar{y}_i +
#' N + \bar{y}^2
#'   }
#' \sqrt{
#' \sum_{i=1}^G n_i \mu_i^2 -
#' 2 \bar{y} \sum_{i=1}^G n_i \mu_i +
#' N + \bar{y}^2
#' }
#' } \right)^2
#' }
#'
#' In this formula, there are \eqn{G} groups (reported observations). (For CvTdb
#' data, a "group" is a specific combination of chemical, species, route,
#' medium, dose, and timepoint.) \eqn{n_i} is the number of subjects for group
#' \eqn{i}. \eqn{\bar{y}_i} is the sample mean for group \eqn{i}. \eqn{s_i} is
#' the sample standard deviation for group \eqn{i}.\eqn{\mu_i} is the
#' model-predicted value for group \eqn{i}.

#' \eqn{\bar{y}} is the grand mean of observations:
#'
#' \deqn{ \bar{y} = \frac{ \sum_{i=1}^G n_i \bar{y}_i }{\sum_{i=1}^G n_i} }
#'
#' \eqn{\bar{\mu}} is the grand mean of predictions:
#'
#' \deqn{ \bar{\mu} = \frac{ \sum_{i=1}^G n_i \mu_i }{\sum_{i=1}^G n_i} }
#'
#' \eqn{N} is the grand total of subjects:
#'
#' \deqn{N = \sum_{i=1}^G n_i}
#'
#' For the non-summary case (\eqn{N} single-subject observations, with all
#' \eqn{n_i = 1}, \eqn{s_i = 0}, and \eqn{\bar{y}_i = y_i}), this formula
#' reduces to the familiar formula
#'
#' \deqn{ r^2 = \left( \frac{\sum_{i=1}^N (y_i - \bar{y}) (\mu_i - \bar{\mu})}
#' {\sqrt{ \sum_{i=1}^N (y_i - \bar{y})^2 }
#' \sqrt{ \sum_{i=1}^N (\mu_i - \bar{\mu})^2 }
#'  } \right)^2
#'  }
#'
#' # Left-censored data
#'
#' If the observed value is censored, and the predicted value is less than the
#' reported LOQ, then the observed value is (temporarily) set equal to the
#' predicted value, for an effective error of zero.
#'
#' If the observed value is censored, and the predicted value is greater than
#' the reported LOQ, the the observed value is (temporarily) set equal to the
#' reported LOQ, for an effective error of (LOQ - predicted).
#'
#' # Log-10 transformation
#'
#' If `log10 %in% TRUE`, then both the observed and predicted values will be
#'  log10-transformed before calculating the RMSE. In the case where
#' observed values are reported in summary format, each sample mean and sample
#' SD (reported on the natural scale, i.e. the mean and SD of natural-scale
#' individual observations) are used to produce an estimate of the log10-scale
#' sample mean and sample SD (i.e., the mean and SD of log10-transformed
#' individual observations), using [convert_summary_to_log10()].
#'
#' The formulas are as follows. Again, \eqn{\bar{y}_i} is the sample mean for
#' group \eqn{i}. \eqn{s_i} is the sample standard deviation for group \eqn{i}.
#'
#' \deqn{\textrm{log10-scale sample mean}_i = \log_{10}
#' \left(\frac{\bar{y}_i^2}{\sqrt{\bar{y}_i^2 + s_i^2}} \right)}
#'
#' \deqn{\textrm{log10-scale sample SD}_i = \sqrt{\log_{10} \left(1 +
#' \frac{s_i^2}{\bar{y}_i^2} \right)}}
#'
#' @param pred Numeric vector: Model-predicted value corresponding to each
#'   observed value. Even if `log10_trans %in% TRUE`, these should *not* be
#'   log-transformed.  (If `log10_trans %in% TRUE`, they will be
#'   log10-transformed internally to this function before calculation.)
#' @param obs Numeric vector: Observed sample means for summary data, or
#'   observed values for non-summary data. Censored observations should *not* be
#'   NA; they should be substituted with the LOQ. Even if `log10_trans %in%
#'   TRUE`, these should *not* be log10-transformed. (If `log10_trans %in% TRUE`,
#'   they will be transformed to log10-scale means internally to this function
#'   before calculation.)
#' @param obs_sd Numeric vector: Observed sample SDs for summary data. For
#'   non-summary data (individual-subject observations), the corresponding
#'   element of `obs_sd` should be set to 0. Even if `log10_trans %in% TRUE`,
#'   these should *not* be log10-transformed. (If `log10_trans %in% TRUE`, they
#'   will be transformed to log10-scale standard deviations internally to this
#'   function before calculation.)
#' @param n_subj Numeric vector: Observed sample number of subjects for summary
#'   data. For non-summary data (individual-subject observations), `group_n`
#'   should be set to 1.
#' @param detect Logical: Whether each
#' @param log10_trans Logical. FALSE (default) means that R-squared is computed for
#'   observations vs. predictions. TRUE means that R-squared is computed for
#'   log10(observations) vs. log10(predictions) (see Details).
#' @return A numeric scalar: the R-squared value for observations vs.
#'   predictions.
#' @author Caroline Ring
calc_rsq <- function(pred,
                     obs,
                     obs_sd,
                     n_subj,
                     detect,
                     log10_trans = FALSE) {
  # If both obs and pred are below LOQ, set pred = LOQ.
  # This will effectively make error zero in these cases.
  pred <- ifelse(detect %in% FALSE & pred <= obs, obs, pred)

  # Convert to log10-scale if necessary
  if (log10_trans %in% TRUE) {
    tmplist <- convert_summary_to_log10(sample_mean = obs,
                                        sample_SD = obs_sd)
    obs <- tmplist$log10mean
    obs_sd <- tmplist$log10SD
    pred <- log10(pred)
  }

  # grand mean of predictions
  pred_bar <- sum(n_subj * pred) / sum(n_subj)
  # grand mean of observations
  x_bar <- sum(n_subj * obs) / sum(n_subj)
  # Calculate numerator of Pearson correlation coefficient
  r_num <- (sum(pred * n_subj * obs)
            - pred_bar * sum(n_subj * obs)
            - x_bar * sum(n_subj * pred)
            + x_bar * pred_bar * sum(n_subj)
  )

  # Calculate first denominator term of Pearson correlation coefficient
  r_denom_1 <- sqrt(sum((n_subj - 1) * obs_sd^2 + n_subj * obs^2)
                    - 2 * x_bar * sum(n_subj * obs)
                    + sum(n_subj) * x_bar^2)

  # Calculate second denominator term of Pearson correlation coefficient
  r_denom_2 <- sqrt(pmax(0, sum(n_subj * pred^2)
                         - 2 * pred_bar * sum(n_subj * pred)
                         + sum(n_subj) * pred_bar^2))

  # Put the pieces together to get Pearson correlation coefficient
  if (isTRUE(r_denom_1 > 0) && isTRUE(r_denom_2 > 0)) {
    r <- r_num / (r_denom_1 * r_denom_2)
  } else {
    r <- NA_real_
  }

  # Square it
  rsq <- r^2

  return(rsq)
}
