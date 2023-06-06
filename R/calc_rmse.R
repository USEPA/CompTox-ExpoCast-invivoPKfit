#' Calculate RMSE (root mean squared error)
#'
#' Calculate RMSE when observed data may be left-censored (non-detect) or may be
#' reported in summary form (as sample mean, sample standard deviation, and
#' sample number of subjects). Additionally, handle the situation when observed
#' data and predictions need to be log-transformed before RMSE is calculated.
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
#' # Left-censored data
#'
#' If the observed value is censored, and the predicted value is less than the
#' reported LOQ, then the observed value is (temporarily) set equal to the
#' predicted value, for an effective error of zero.
#'
#' If the observed value is censored, and the predicted value is greater than
#' the reported LOQ, the the observed value is set equal to the reported LOQ.
#'
#' # Log transformation
#'
#' If `log %in% TRUE`, then both the observed and predicted values will be
#' (natural) log-transformed before calculating the RMSE. In the case where
#' observed values are reported in summary format, each sample mean and sample
#' SD (reported on the natural scale, i.e. the mean and SD of natural-scale
#' individual observations) are used to produce an estimate of the log-scale
#' sample mean and sample SD (i.e., the mean and SD of log-transformed
#' individual observations), using [convert_summary_to_log10()].
#'
#' The formulas are as follows. Again, \eqn{\bar{y}_i} is the sample mean for
#' group \eqn{i}. \eqn{s_i} is the sample standard deviation for group \eqn{i}.
#'
#' \deqn{\textrm{log-scale sample mean}_i = \log
#' \left(\frac{\bar{y}_i^2}{\sqrt{\bar{y}_i^2 + s_i^2}} \right)}
#'
#' \deqn{\textrm{log-scale sample SD}_i = \sqrt{\log \left(1 +
#' \frac{s_i^2}{\bar{y}_i^2} \right)}}
#'
#' @param group_mean Numeric vector: Observed sample means for summary data, or
#'   observed values for non-summary data. Censored observations should *not* be
#'   NA; they should be substituted with some value at or below the
#'   corresponding LOQ (e.g. LOQ or LOQ/2). Even if `log %in% TRUE`, these
#'   should *not* be log-transformed. (If `log %in% TRUE`, they will be
#'   transformed to log-scale means internally to this function before
#'   calculation.)
#' @param group_sd Numeric vector: Observed sample SDs for summary data. For
#'   non-summary data (individual-subject observations), the corresponding
#'   element of `group_sd` should be set to 0. Even if `log %in% TRUE`, these
#'   should *not* be log-transformed. (If `log %in% TRUE`, they will be
#'   transformed to log-scale standard deviations internally to this function
#'   before calculation.)
#' @param group_n Numeric vector: Observed sample number of subjects for summary
#'   data. For non-summary data (individual-subject observations), `group_n`
#'   should be set to 1.
#' @param pred Numeric vector: Model-predicted value corresponding to each
#'   observed value. Even if `log %in% TRUE`, these should *not* be
#'   log-transformed. (If `log %in% TRUE`, they will be log-transformed
#'   internally to this function before calculation.)
#' @param group_LOQ Numeric vector: Reported limit of quantification (LOQ) (i.e,
#'   left-censoring limit) for each observation. May contain `NA_real_`s for
#'   non-censored observations. Even if `log %in% TRUE`, these should *not* be
#'   log-transformed. (If `log %in% TRUE`, log transformation will be handled
#'   internally to this function before calculation.)
#' @param log Logical. FALSE (default) means that RMSE is computed for
#'   natural-scale observations and predictions -- effectively, `sqrt(mean(
#'   (observed - pred)^2 ))`. TRUE means that observations and predictions will
#'   be log-transformed before RMSE is calculated (see Details)
#'   -- effectively `sqrt(mean( ( log(observed) - log(pred) )^2 ))`.
#' @return A numeric scalar: the root mean squared error (RMSE) for this set of
#'   observations and predictions.
#' @author Caroline Ring
calc_rmse <- function(pred,
                      obs,
                      obs_sd,
                      n_subj,
                      detect,
                      log10_trans = FALSE){
  #If both obs and pred are below LOQ, set pred to obs
  #This will effectively make error zero in these cases.
  pred <- ifelse(detect %in% FALSE &
                   (pred <= obs) %in% TRUE,
                 obs,
                 pred)

  #Convert to log10-scale if necessary
  if(log10_trans %in% TRUE){
    tmplist <- convert_summary_to_log10(sample_mean = obs,
                                      sample_SD = obs_sd)
    obs <- tmplist$log10mean
    obs_sd <- tmplist$log10SD
    pred <- log10(pred)
  }

#Calculate MSE.
  #If log = TRUE, this is mean(( log(observed) - log(pred) )^2)
  #If log = FALSE, this is mean((observed - pred)^2)
  mse <- (1/sum(n_subj)) *
    sum(
    (n_subj - 1) * obs_sd^2 +
      n_subj * obs^2  +
      -2 * pred * n_subj * obs  +
      n_subj * pred^2
  )

  #RMSE is just sqrt of MSE
  rmse <- sqrt(mse)

  return(rmse)
}

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
#' # Log transformation
#'
#' If `log %in% TRUE`, then both the observed and predicted values will be
#' (natural) log-transformed before calculating the RMSE. In the case where
#' observed values are reported in summary format, each sample mean and sample
#' SD (reported on the natural scale, i.e. the mean and SD of natural-scale
#' individual observations) are used to produce an estimate of the log-scale
#' sample mean and sample SD (i.e., the mean and SD of log-transformed
#' individual observations), using [convert_summary_to_log10()].
#'
#' The formulas are as follows. Again, \eqn{\bar{y}_i} is the sample mean for
#' group \eqn{i}. \eqn{s_i} is the sample standard deviation for group \eqn{i}.
#'
#' \deqn{\textrm{log-scale sample mean}_i = \log
#' \left(\frac{\bar{y}_i^2}{\sqrt{\bar{y}_i^2 + s_i^2}} \right)}
#'
#' \deqn{\textrm{log-scale sample SD}_i = \sqrt{\log \left(1 +
#' \frac{s_i^2}{\bar{y}_i^2} \right)}}
#'
#' @param obs Numeric vector: Observed sample means for summary data, or
#'   observed values for non-summary data. Censored observations should *not* be
#'   NA; they should be substituted with the LOQ. Even if `log10_trans %in%
#'   TRUE`, these should *not* be log-transformed. (If `log10_trans %in% TRUE`,
#'   they will be transformed to log-scale means internally to this function
#'   before calculation.)
#' @param obs_sd Numeric vector: Observed sample SDs for summary data. For
#'   non-summary data (individual-subject observations), the corresponding
#'   element of `obs_sd` should be set to 0. Even if `log10_trans %in% TRUE`,
#'   these should *not* be log-transformed. (If `log10_trans %in% TRUE`, they
#'   will be transformed to log-scale standard deviations internally to this
#'   function before calculation.)
#' @param n_subj Numeric vector: Observed sample number of subjects for summary
#'   data. For non-summary data (individual-subject observations), `group_n`
#'   should be set to 1.
#' @param detect Logical: Whether each
#' @param pred Numeric vector: Model-predicted value corresponding to each
#'   observed value. Even if `log10_trans %in% TRUE`, these should *not* be
#'   log-transformed.  (If `log10_trans %in% TRUE`, they will be
#'   log10-transformed internally to this function before calculation.)
#' @param log Logical. FALSE (default) means that R-squared is computed for
#'   observations vs. predictions. TRUE means that R-squared is computed for
#'   log(observations) vs. log(predictions) (see Details).
#' @return A numeric scalar: the R-squared value for observations vs.
#'   predictions.
#' @author Caroline Ring
calc_rsq <- function(obs,
                     obs_sd,
                     n_subj,
                     detect,
                     pred,
                     log10_trans = FALSE){
  #If both obs and pred are below LOQ, set pred = LOQ.
  #This will effectively make error zero in these cases.
  pred <- ifelse(detect %in% FALSE &
                  pred <= obs,
                obs,
                pred)

  #Convert to log10-scale if necessary
  if(log10_trans %in% TRUE){
    tmplist <- convert_summary_to_log10(sample_mean = obs,
                                      sample_SD = obs_sd)
    obs <- tmplist$log10mean
    obs_sd <- tmplist$log10SD
    pred <- log10(pred)
  }

  #grand mean of predictions
  pred_bar <- sum(n_subj*pred)/sum(n_subj)
  #grand mean of observations
  x_bar <- sum(n_subj*obs)/sum(n_subj)
#Calculate numerator of Pearson correlation coefficient
  r_num <- sum(pred * n_subj * obs) -
    pred_bar * sum(n_subj * obs) -
    x_bar * sum(n_subj * pred) +
    x_bar * pred_bar * sum(n_subj)
#Calculate first denominator term of Pearson correlation coefficient
  r_denom_1 <- sqrt(sum(
    (n_subj - 1) * obs_sd^2 +
      n_subj * obs^2
  ) -
    2 * x_bar * sum(n_subj * obs) +
    sum(n_subj) * x_bar^2)
  #Calculate second denominator term of Pearson correlation coefficient
  r_denom_2 <- sqrt(sum(n_subj * pred^2) -
    2 * pred_bar * sum(n_subj * pred) +
    sum(n_subj) * pred_bar^2)

  #Put the pieces together to get Pearson correlation coefficient
  if((r_denom_1 > 0) %in% TRUE  &
    ( r_denom_2 > 0) %in% TRUE){
    r <- r_num /(r_denom_1 * r_denom_2)
  }else{
    r <- NA_real_
  }

#Square it
rsq <- r^2

return(rsq)

}

#' Combined standard deviation
#'
#' Given mean, standard deviation, and N for some set of groups, calculate the
#' combined standard deviation. Note that the groups may not overlap.
#'
#' @param group_mean Numeric vector: Observed sample means for summary data, or
#'   observed values for non-summary data. Censored observations should *not* be
#'   NA; they should be substituted with some value at or below the
#'   corresponding LOQ (e.g. LOQ or LOQ/2). Even if `log %in% TRUE`, these
#'   should *not* be log-transformed.
#' @param group_sd Numeric vector: Observed sample SDs for summary data. For
#'   non-summary data (individual-subject observations), the corresponding
#'   element of `group_sd` should be set to 0. Even if `log %in% TRUE`, these
#'   should *not* be log-transformed.
#' @param group_n Numeric vector: Observed sample number of subjects for summary
#'   data. For non-summary data (individual-subject observations), `group_n`
#'   should be set to 1.
#' @param unbiased Logical. If TRUE (the default), then `group_sd` is assumed to
#'   be the unbiased estimator of population standard deviation (i.e. calculated
#'   using `n-1` in the denominator -- the way that `stats::sd()` calculates
#'   it), and the returned combined SD is also the unbiased estimator of the
#'   combined population SD. If FALSE, then `group_sd` is assumed to be the
#'   biased estimator (using `n` in the denominator), and the returned value is
#'   also the biased estimator of the combined population SD.
#' @param na.rm Logical. If TRUE (default), then any groups where mean, SD, *or*
#'   N were NA will be dropped. If FALSE, they will be retained (and the result
#'   will be NA).
#' @return Numeric: the standard deviation of the combined population (i.e. if
#'   all the groups were concatenated into one large group).
#' @author Caroline Ring
combined_sd <- function(group_mean,
                        group_sd,
                        group_n,
                        unbiased = TRUE,
                        na.rm = TRUE,
                        log = FALSE){

  x_len <- c("group_mean" = length(group_mean),
             "group_sd" = length(group_sd),
             "group_n" = length(group_n))

  if(any(x_len %in% 0)){
    stop(paste0("invivopkfit::combined_sd(): ",
                "the following arguments have zero length: ",
                paste(names(x_len)[x_len %in% 0],
                      collapse = ", ")
    ))
  }

  max_len <- max(x_len)
  which_max_len <- which.max(x_len)

  bad_len <- (x_len < max_len) & (x_len != 1)


  if(any(bad_len)){
    warning(paste("invivopkfit::combined_sd():",
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
    ))
  }

  #repeat to match longest
  for (i in seq_along(x_len)){
    assign(names(x_len)[i],
           rep( #repeat the current value of each item to match the length
             get(names(x_len)[i]), #get the current value of each item
             length.out = max_len)
    )
  }

  #remove NAs if so specified
  if(na.rm %in% TRUE){
    which_na <- is.na(group_mean) | is.na(group_sd) | is.na(group_n)
    group_mean <- group_mean[!which_na]
    group_sd <- group_sd[!which_na]
    group_n <- group_n[!which_na]
  }

  grand_mean <- sum(group_n*group_mean)/sum(group_n)

  #if all N = 1, then just take regular standard deviation
  if(all(group_n %in% 1)){
    grand_sd <- sd(group_mean)

    grand_n <- sum(group_n)
    if(unbiased %in% FALSE){
      #convert unbiased SD to biased SD
      grand_sd <- grand_sd * sqrt((grand_n-1)/grand_n)
    }
  }else{ #if not all N = 1
    if(unbiased %in% TRUE){
      #convert unbiased group SDs to biased group SDs
      group_sd <- group_sd *
        sqrt((group_n[which_na]-1)/group_n)
    }

    grand_var <- (sum(group_n*group_sd^2) +
                    sum(group_n*(group_mean - grand_mean)^2))/
      (sum(group_n))

    grand_sd <- sqrt(grand_var) #biased


    if(unbiased %in% TRUE){
      grand_n <- sum(group_n)
      #convert biased grand SD to unbiased grand SD
      grand_sd <- grand_sd * sqrt(grand_n/(grand_n - 1))
    }
  }

  if(log %in% TRUE){
    #convert to log-scale combined SD
    grand_sd <- sqrt(log(1 + grand_sd^2 / grand_mean^2))
  }

  return(grand_sd)
}

#' Combined mean
#'
#' Given mean, standard deviation, and N for some set of groups, calculate the
#' combined mean. Note that the groups may not overlap.
#'
#' @param group_mean Numeric vector: Observed sample means for summary data, or
#'   observed values for non-summary data. Censored observations should *not* be
#'   NA; they should be substituted with some value at or below the
#'   corresponding LOQ (e.g. LOQ or LOQ/2). Even if `log %in% TRUE`, these
#'   should *not* be log-transformed.
#' @param group_n Numeric vector: Observed sample number of subjects for summary
#'   data. For non-summary data (individual-subject observations), `group_n`
#'   should be set to 1.
#' @param na.rm Logical. If TRUE (default), then any groups where mean *or*
#'   N were NA will be dropped. If FALSE, they will be retained (and the result
#'   will be NA).
#' @return Numeric: the mean of the combined population (i.e. if
#'   all the groups were concatenated into one large group).
#' @author Caroline Ring
combined_mean <- function(group_mean,
                        group_n,
                        na.rm = TRUE,
                        log = FALSE){

  x_len <- c("group_mean" = length(group_mean),
             "group_n" = length(group_n))

  if(any(x_len %in% 0)){
    stop(paste0("invivopkfit::combined_mean(): ",
                "the following arguments have zero length: ",
                paste(names(x_len)[x_len %in% 0],
                      collapse = ", ")
    ))
  }

  max_len <- max(x_len)
  which_max_len <- which.max(x_len)

  bad_len <- (x_len < max_len) & (x_len != 1)


  if(any(bad_len)){
    warning(paste("invivopkfit::combined_sd():",
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
    ))
  }

  #repeat to match longest
  for (i in seq_along(x_len)){
    assign(names(x_len)[i],
           rep( #repeat the current value of each item to match the length
             get(names(x_len)[i]), #get the current value of each item
             length.out = max_len)
    )
  }

  #remove NAs if so specified
  if(na.rm %in% TRUE){
    which_na <- is.na(group_mean) | is.na(group_sd) | is.na(group_n)
    group_mean <- group_mean[!which_na]
    group_n <- group_n[!which_na]
  }

  grand_mean <- sum(group_n*group_mean)/sum(group_n)

  return(grand_mean)
}

#' Breusch-Pagan test of heteroscedasticity
#' Evaluate heteroscedasticity using Breusch-Pagan test
#' The Breusch-Pagan test performs a linear regression of squared residuals on one or more explanatory variables. The test statistic is formed by the R-squared value of this regression multiplied by the number of residuals. This test statistic obeys a chi-squared distribution with S-1 degrees of freedom, where S is the number of explanatory terms (plus an intercept). The p-value of the test is given accordingly.
#'
#' @param varformula A formula for the variance linear model. The left-hand side of this formula must be the squared residuals.
#' @param data A `data.frame` containing all variables referenced in `varformula`.
#' @return A list with the following elements:
#' - `statistic` : The test statistic
#' - `p.value`: The p-value of the test
#' - `r.squared`: The R-squared of the linear regression described in `varformula`
#' - `n`: The number of observations
#' - `df`: The chi-squared degrees of freedom
#' @author Caroline Ring
#' @references https://rpubs.com/cyobero/187387

bp_test <- function(group_mean,
                    group_sd,
                    group_n,
                    group_LOQ,
                    pred,
                    log = FALSE){
  n_obs <- length(group_mean)
  #check if there are more than 1 unique predicted value to 6 places
  #i.e., whether this is a flat model
  #if it's flat, don't bother calculating,
  #as we know rsq of resid^2 ~ pred will be 0.
  do_calc <- length(unique(signif(pred, 6)))>1

  if(do_calc){
  #If both obs and pred are below LOQ, set obs to pred.
  #This will effectively make error zero in these cases.
  if(any((group_mean <= group_LOQ &
          pred <= group_LOQ) %in% TRUE)){
    group_mean[(group_mean <= group_LOQ &
                  pred <= group_LOQ) %in% TRUE] <- pred[(group_mean <= group_LOQ &
                                                           pred <= group_LOQ) %in% TRUE]
  }

  #If obs is below LOQ but pred is not, set obs to LOQ.
  if(any((group_mean <= group_LOQ &
          !(pred <= group_LOQ)) %in% TRUE)){
    group_mean[(group_mean <= group_LOQ &
                  !(pred <= group_LOQ)) %in% TRUE] <- group_LOQ[(group_mean <= group_LOQ &
                                                                   !(pred <= group_LOQ)) %in% TRUE]
  }

  #Convert to log-scale if necessary
  if(log %in% TRUE){
    tmplist <- convert_summary_to_log10(sample_mean = group_mean,
                                      sample_SD = group_sd)
    group_mean <- tmplist$logmean
    group_sd <- tmplist$logSD
    pred <- log(pred)
  }

  #we need to fit b0 and b1 for a linear model of squared residuals vs. pred

    #if all N_Subjects = 1, we can use lm(); otherwise we need to use optimx()
    #so that we can use summary log-likelihood
    if(all(group_n %in% 1)){
      resid2 <- (group_mean - pred)^2
      bp_lm <- lm(resid2 ~ pred)
      bp_rsq <- summary(bp_lm)$r.squared
    }else{
  opt_par <- tryCatch({
    tmp <- optimx::optimx(par = c("b0" = 0,
                                "b1" = 0,
                                "sigma" = 1),
                        fn = bp_loglike_summary,
                        method = "bobyqa",
                        lower = c("b0" = -1000,
                                  "b1" = -1000,
                                  "sigma" = 1e-8),
                        upper = c("b0" = 1e4,
                                  "b1" = 1e4,
                                  "sigma" = 1e4),
                        Ni = group_n,
                        ybari = group_mean,
                        si = group_sd,
                        mui = pred)
  stats::coefficients(tmp)[1,]
  },
  error = function(err) return(c("b0" = NA_real_, "b1" = NA_real_, "sigma" = NA_real_))
  )
  b0 <- opt_par["b0"]
  b1 <- opt_par["b1"]


  bp_rsq <- tryCatch(
    1 - bp_ssr(b0 = b0,
                      b1 = b1,
                      Ni = group_n,
                      ybari = group_mean,
                      si = group_sd,
                      mui = pred)/bp_sst(Ni = group_n,
                                         ybari = group_mean,
                                         si = group_sd,
                                         mui = pred),
    error = function(err) return(NA_real_))
    }


  qval <- bp_rsq * n_obs
  pval <- pchisq(q = qval,
                 df = 1,
                 lower.tail = FALSE)
  }else{
    bp_rsq <- 0
    qval <- 0
    pval <- 1
  }
  return(list("statistic" = qval,
              "p.value" = pval,
              "r.squared" = bp_rsq,
              "n" = n_obs,
              "df" = 1))
}

bp_loglike_summary <- function(par, Ni, ybari, si, mui, negative = TRUE){
  b0 = par["b0"]
  b1 = par["b1"]
  sigma = par["sigma"]

  N <- sum(Ni)
  LL <- N * log(1/(sigma*sqrt(1*pi))) -
    (1/(2*sigma^2)) * bp_ssr(b0 = b0,
                             b1 = b1,
                             Ni = Ni,
                             ybari = ybari,
                             si = si,
                             mui = mui)

  if(negative %in% TRUE){
    LL <- -1 * LL
  }

  return(LL)
}

bp_ssr <- function(b0, b1,
                   Ni,
                   ybari,
                   si,
                   mui){
  N <- sum(Ni)
  sum(3*N*si^4 +
        2*si^2*(b0 + b1*mui -
                  3*(ybari - mui)^2) +
        (b0^2 + 2*b0*(b1*mui -
                        si^2 -
                        (ybari - mui)^2) +
           b1^2*mui^2 +
           2*b1*(-si^2 -
                   (ybari - mui)^2)*mui +
           6*si^2*mui^2 +
           ybari^4 -
           4*ybari^3*mui +
           ybari^2*(6*si^2 + 6*mui^2) +
           ybari*(-12*si^2*mui - 4*mui^3) +
           mui^4)*Ni)
}

bp_sst <- function(Ni,
                   ybari,
                   si,
                   mui){
  N <- sum(Ni)
  #eps2bar = grand mean of squared residuals
  eps2bar <- (1/N)*sum((Ni-1) * si^2 +
                         Ni * ybari^2 -
                         2 * mui * Ni * ybari +
                         Ni * mui^2)

  sum(3*N*si^4 +
        eps2bar^2*Ni +
        4*eps2bar*ybari*Ni*mui -
        2*eps2bar*(si^2*(Ni - 1) +
                     ybari^2*Ni) -
        2*eps2bar*Ni*mui^2 +
        6*si^2*ybari^2*Ni -
        6*si^2*ybari^2 +
        ybari^4*Ni -
        4*ybari*(3*si^2*Ni -
                   3*si^2 +
                   ybari^2*Ni)*mui -
        4*ybari*Ni*mui^3 +
        6*(si^2*(Ni - 1) +
             ybari^2*Ni)*mui^2 +
        Ni*mui^4)
}

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
#'
#' @param sample_mean Numeric: one or more sample means
#' @param sample_SD Numeric: one or more sample SDs
#' @return A list with two named elements: "log10mean" and "log10SD", the log10-scale
#'   sample means and log10-scale sample SDs, respectively.
#' @export
#' @author Caroline Ring
convert_summary_to_log10 <- function(sample_mean, sample_SD){
  log10mean <- log10(sample_mean^2 / sqrt(sample_SD^2 + sample_mean^2))
  log10SD <- sqrt(log10(1 + (sample_SD^2 / sample_mean^2)))
  return(list("log10mean" = log10mean,
              "log10SD" = log10SD))
}
