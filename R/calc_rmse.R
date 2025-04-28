#' Calculate RMSE (root mean squared error)
#'
#' Calculate RMSE when observed data may be left-censored (non-detect) or may be
#' reported in summary form (as sample mean, sample standard deviation, and
#' sample number of subjects). Additionally, handle the situation when observed
#' data and predictions need to be log10-transformed before RMSE is calculated.
#'
#' RMSE is calculated using the following formula, to properly handle summary
#' data:
#'
#' \deqn{
#' \sqrt{
#' \frac{1}{N}
#'  \sum_{i=1}^G \left( (n_i - 1) s_i^2 +
#'   n_i \bar{y}_i^2 - 2 n_i \bar{y}_i \mu_i + \mu_i^2 \right)
#'    }
#' }
#'
#' In this formula, there are \eqn{G} groups. (For CvTdb data, a "group" is a
#' specific combination of chemical, species, route, medium, dose, and
#' timepoint.) \eqn{n_i} is the number of subjects for group \eqn{i}.
#' \eqn{\bar{y}_i} is the sample mean for group \eqn{i}. \eqn{s_i} is the sample
#' standard deviation for group \eqn{i}.\eqn{\mu_i} is the model-predicted value
#' for group \eqn{i}.
#'
#' \eqn{N} is the grand total of subjects:
#'
#' \deqn{N = \sum_{i=1}^G n_i}
#'
#' For the non-summary case (\eqn{N} single-subject observations, with all
#' \eqn{n_i = 1}, \eqn{s_i = 0}, and \eqn{\bar{y}_i = y_i}), this formula
#' reduces to the familiar RMSE formula
#'
#' \deqn{\sqrt{\frac{1}{N} \sum_{i=1}^N (y_i - \mu_i)^2}}
#'
#' @section Left-censored data:
#'
#' If the observed value is censored, and the predicted value is less than the
#' reported LOQ, then the observed value is (temporarily) set equal to the
#' predicted value, for an effective error of zero.
#'
#' If the observed value is censored, and the predicted value is greater than
#' the reported LOQ, the the observed value is set equal to the reported LOQ.
#'
#' @section Log10 transformation:
#'
#' If `log10_trans %in% TRUE`, then both the observed and predicted values will be
#' log10-transformed before calculating the RMSE. In the case where
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
#'   log-transformed. (If `log10_trans %in% TRUE`, they will be log-transformed
#'   internally to this function before calculation.)
#' @param obs Numeric vector: Observed sample means for each observation if
#'   summary data, or observed values for each observation if non-summary data.
#'   Censored observations should *not* be NA; they should be substituted with
#'   the LOQ. Even if `log10_trans %in% TRUE`, these should *not* be
#'   log-transformed. (If `log10_trans %in% TRUE`, they will be transformed to
#'   log-scale means internally to this function before calculation.)
#' @param obs_sd Numeric vector: Observed sample SDs for each observation, if
#'   summary data. For non-summary data (individual-subject observations), the
#'   corresponding element of `group_sd` should be set to 0. Even if
#'   `log10_trans %in% TRUE`, these should *not* be log-transformed. (If
#'   `log10_trans %in% TRUE`, they will be transformed to log-scale standard
#'   deviations internally to this function before calculation.)
#' @param n_subj Numeric vector: Observed sample number of subjects for each
#'   observation. For non-summary data (individual-subject observations),
#'   `n_subj` should be set to 1.
#' @param detect Logical vector: `TRUE` for each observation that was detected
#'   (above LOQ); `FALSE` for each observation that was non-detect (below LOQ).
#' @param log10_trans Logical. FALSE (default) means that RMSE is computed for
#'   natural-scale observations and predictions -- effectively, `sqrt(mean(
#'   (observed - pred)^2 ))`. TRUE means that observations and predictions will
#'   be log10-transformed before RMSE is calculated (see Details)
#'   -- effectively `sqrt(mean( ( log(observed) - log(pred) )^2 ))`.
#' @return A numeric scalar: the root mean squared error (RMSE) for this set of
#'   observations and predictions.
#' @author Caroline Ring
calc_rmse <- function(pred,
                      obs,
                      obs_sd,
                      n_subj,
                      detect,
                      log10_trans = FALSE) {
  # If both obs and pred are below LOQ, set pred to obs
  # This will effectively make error zero in these cases.
  pred <- ifelse(detect %in% FALSE & (pred <= obs) %in% TRUE, obs, pred)

  # Convert to log10-scale if necessary
  if (isTRUE(log10_trans)) {
    tmplist <- convert_summary_to_log10(sample_mean = obs,
                                      sample_SD = obs_sd)
    obs <- tmplist$log10mean
    obs_sd <- tmplist$log10SD
    pred <- log10(pred)
  }

# Calculate MSE.
  # If log = TRUE, this is mean(( log(observed) - log(pred) )^2)
  # If log = FALSE, this is mean((observed - pred)^2)
  mse <- (1 / sum(n_subj)) * sum((n_subj - 1) * obs_sd^2
                                 + n_subj * obs^2
                                 - 2 * pred * n_subj * obs
                                 + n_subj * pred^2)

  # RMSE is just sqrt of MSE
  rmse <- sqrt(mse)

  return(rmse)
}
