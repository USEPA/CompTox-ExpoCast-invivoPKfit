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

