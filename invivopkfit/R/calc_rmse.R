calc_rmse <- function(group_mean,
                      group_sd,
                      group_N,
                      group_pred,
                      group_LOQ,
                      log = FALSE){
  #If both obs and pred are below LOQ, set obs to pred.
  #This will effectively make error zero in these cases.
  group_mean[group_mean <= group_LOQ &
               group_pred <= group_LOQ] <- group_pred[group_mean < group_LOQ &
                                                       group_pred < group_LOQ]

  #If obs is below LOQ but pred is not, set obs to LOQ.
  group_mean[group_mean <= group_LOQ &
             !(group_pred <= group_LOQ)] <- group_LOQ[group_mean <= group_LOQ &
                                                       !(group_pred <= group_LOQ)]

  #Convert to log-scale if necessary
  if(log %in% TRUE){
    tmplist <- convert_summary_to_log(sample_mean = group_mean,
                                      sample_SD = group_sd)
    group_mean <- tmplist$logmean
    group_sd <- tmplist$logSD
  }

#Calculate MSE.
  mse <- (1/sum(group_N)) * sum(
    (group_N - 1) * group_sd^2 +
      group_N * group_mean^2 +
      -2 * group_N * group_mean * group_pred +
      group_pred^2
  )

  rmse <- sqrt(mse)

  return(rmse)
}

calc_rsq <- function(group_mean,
                     group_sd,
                     group_N,
                     group_pred,
                     log = FALSE){
  #If both obs and pred are below LOQ, set obs to pred.
  #This will effectively make error zero in these cases.
  group_mean[group_mean <= group_LOQ &
               group_pred <= group_LOQ] <- group_pred[group_mean < group_LOQ &
                                                       group_pred < group_LOQ]

  #If obs is below LOQ but pred is not, set obs to LOQ.
  group_mean[group_mean <= group_LOQ &
               !(group_pred <= group_LOQ)] <- group_LOQ[group_mean < group_LOQ &
                                                         !(group_pred < group_LOQ)]

  #Convert to log-scale if necessary
  if(log %in% TRUE){
    tmplist <- convert_summary_to_log(sample_mean = group_mean,
                                      sample_SD = group_sd)
    group_mean <- tmplist$logmean
    group_sd <- tmplist$logSD
  }

  pred_bar <- weighted.mean(group_pred, group_N)
  x_bar <- weighted.mean(group_mean, group_N)
#Calculate numerator of Pearson correlation coefficient
  r_num <- sum(group_pred * group_N * group_mean) -
    pred_bar * sum(group_N * group_mean) -
    x_bar * sum(group_N * group_pred) +
    x_bar * pred_bar * sum(group_N)
#Calculate first denominator term of Pearson correlation coefficient
  r_denom_1 <- sqrt(sum(
    (group_N - 1) * group_sd^2 +
      group_N * group_mean^2
  ) -
    2 * x_bar * sum(group_N * group_mean) +
    sum(group_N) * x_bar^2)
  #Calculate second denominator term of Pearson correlation coefficient
  r_denom_2 <- sqrt(sum(group_N * group_pred^2) -
    2 * pred_bar * sum(group_N * group_pred) +
    sum(group_N) * pred_bar^2)

r <- r_num /(r_denom_1 * r_denom_2)

rsq <- r^2

return(rsq)

}

#' Convert sample mean and SD to log-scale
#'
#' @param sample_mean Numeric: one or more sample means
#' @param sample_SD Numeric: one or more sample SDs
#' @return A list with two named elements: "logmean" and "logSD", the log-scale
#'   sample means and log-scale sample SDs, respectively.
#' @author Caroline Ring
convert_summary_to_log <- function(sample_mean, sample_SD){
  logmean <- log(sample_mean^2 / sqrt(sample_SD^2 + sample_mean^2))
  logSD <- sqrt(log(1 + (sample_SD^2 / sample_mean^2)))
  return(list("logmean" = logmean,
              "logSD" = logSD))
}
