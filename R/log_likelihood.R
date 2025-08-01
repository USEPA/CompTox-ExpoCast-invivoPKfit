#' Log-likelihood
#'
#' The log-likelihood function (probability of data given model parameters).
#'
#' The log-likelihood is formulated by assuming that residuals (transformed
#' model-predicted concentrations minus transformed observed concentrations) are
#' independent, and that groups of residuals obey zero-mean normal distributions
#' where each group may have a separate error variance. Error groups are defined
#' as unique combinations of variables in the harmonized data, by a command such
#' as `pk(data = ...) + stat_error_model(error_group = vars(...)`.
#'
#' @section Log-likelihood equations:
#'
#' For chemical-species combination \eqn{i} and study \eqn{j}, define the
#' following quantities.
#'
#' \eqn{y_{ijk}} is the \eqn{k^{th}} observation of concentration (which may be
#' transformed, e.g. dose-normalized and/or log10-transformed), corresponding to
#' dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. Each observation has a
#' corresponding LOQ, \eqn{\textrm{LOQ}_{ijk}}.
#'
#' For multiple-subject observations, \eqn{y_{ijk}} is the \eqn{k^{th}} *sample
#' mean* observed concentration for chemical-species combination \eqn{i} and
#' study \eqn{j}, corresponding to dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. It
#' represents the mean of \eqn{n_{ijk}} individual measurements. It has a
#' corresponding sample standard deviation, \eqn{s_{ijk}}. In the harmonized
#' data,  \eqn{s_{ijk}} is contained in variable `Conc_SD`.
#'
#' \eqn{\bar{\theta}_i} represents the vector of model parameters supplied in
#' argument `params` for chemical-species combination \eqn{i}.
#'
#' \eqn{\mu_{ijk}} is the model-predicted concentration for dose \eqn{d_{ijk}}
#' and time \eqn{t_{ijk}}. If \eqn{f(d, t; \bar{ \theta})} is the model function
#' evaluated at dose \eqn{d} and time \eqn{t}, with parameter vector
#' \eqn{\bar{\theta}}, then
#'
#'  \deqn{\mu_{ijk} = f \left(d_{ijk}, t_{ijk};
#' \bar{\theta}_i \right)}
#'
#' \eqn{\sigma_{ij}^2} is the study- and chemical-specific residual variance. (It
#' is a hyperparameter.)
#'
#' @section Single-subject observations above limit of quantification (detects):
#'
#' This is the normal probability density function evaluated at the observed
#' concentration, as implemented in [stats::dnorm()].
#'
#' \deqn{LL_{ijk} = \log \left( \frac{1}{\sqrt{\sigma_{ij} 2 \pi}} \exp \left[
#' \frac{-1}{2} \left( \frac{y_{ijk} - \mu_{ijk}}{\sigma_{ij}} \right)^2
#' \right]   \right) }
#'
#'
#' @section Single-subject observations below limit of quantification (non-detects):
#'
#' This is the normal cumulative density function evaluated at the LOQ, as
#' implemented in [stats::pnorm()]. It is the total probability of observing a
#' concentration anywhere between 0 and the LOQ.
#'
#' \deqn{LL_{ijk} =  \log \left( \frac{1}{2} \left[ 1 + \textrm{erf} \left(
#' \frac{\textrm{LOQ}_{ijk} - \mu_{ijk}}{\sigma_{ij} \sqrt{2}} \right) \right]  \right) }
#'
#' @section Multiple-subject observations above limit of quantification:
#'
#' This is the joint log-likelihood across the multiple subjects included in one
#' observation, re-expressed in terms of the sample mean, sample SD, and number
#' of subjects for that observation. It is implemented in [dnorm_summary()].
#'
#' \deqn{LL_{ijk} = n_{ijk} * \log \left[
#'  \frac{1}{\sigma_{ij} \sqrt{2 \pi}} +
#'   \frac{-1}{2 \sigma_{ij}^2}  \left(
#'    \left( n_{ijk} - 1 \right)
#'      s_{ijk}^2 + n_{ijk} \left( y_{ijk} - \mu_{ijk} \right)^2
#'       \right)
#'        \right]}
#'
#' @section Multiple-subject observations below limit of quantification:
#'
#' This case is not implemented. If sample mean concentration is reported below
#' LOQ, then it is unclear what individual observed concentrations are
#' represented, and how they were combined to produce the summary data in terms
#' of sample mean and sample SD. Were all individual observations below LOQ? Or
#' were below-LOQ observations replaced with 0, LOQ/2, etc. before sample mean
#' and sample SD were computed? If the sample mean is reported below LOQ, what
#' LOQ is reported? Did individual observations all have the same LOQ, or is an
#' average or median LOQ being used?  It is impossible to formulate the
#' log-likelihood without knowing the answers to these questions. Therefore,
#' multiple-subject observations below LOQ are excluded from analysis (they are
#' marked as excluded in [preprocess_data()]).
#'
#' @section Joint log-likelihood for a chemical and species:
#'
#' The joint log-likelihood for a chemical and species is simply the sum of
#' log-likelihoods across observations.
#'
#' \deqn{LL_{i} = \sum_{j=1}^{J_i} \sum_{k=1}^{K_{ij}} LL_{ijk}}
#'
#' This is the overall probability of the observed data, given the model and
#' parameters.
#'
#' @param par A named list of parameters and their values that are being
#'  optimized.
#' @param const_params Optional: A named list of parameters and their values that
#'  are being held constant.
#' @param data A `data.frame` of data with harmonized variable names. Required:
#'  `Time_trans`, `Dose`, `Conc`, `Detect`, `N_Subjects`, `Conc_SD`. `Conc` and
#'  `Conc_SD` will be transformed according to `dose_norm` and `log10_trans`.
#' @param data_sigma_group A `factor` vector which could be a new variable in
#'  `data`, giving the error group for each row in `data`.
#' @param modelfun Character or function: The name of the function that produces
#'  concentration predictions for the model being evaluated.
#' @param dose_norm Logical: Whether to dose-normalize predicted and observed
#'  concentrations before calculating likelihood.
#' @param log10_trans Logical: Whether to apply a [log10()] transformation to
#'  predicted and observed concentrations before calculating likelihood.
#' @param negative Logical: Whether to return the *negative* log-likelihood
#'  (i.e., the log-likelihood multiplied by negative 1). Default TRUE, to
#'  multiply the log-likelihood by negative 1 before returning it. This option
#'  is useful when treating the log-likelihood as an objective function to be
#'  *minimized* by an optimization algorithm.
#' @param force_finite Logical: Whether to force return of a finite value (e.g.
#'  as required by method `L-BFGS-B` in [optimx::optimx()]). Default FALSE. If
#'  TRUE, then if the log-likelihood works out to be non-finite, then it will be
#'  replaced with `.Machine$double.xmax`.
#' @param includes_preds Logical: whether `data` includes predictions.
#' @param suppress.messages Logical.
#' @return A log-likelihood value for the data given the parameter values in
#'  params
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
log_likelihood <- function(par,
                           const_params = NULL,
                           data = NULL,
                           data_sigma_group = NULL,
                           modelfun = NULL,
                           dose_norm = FALSE,
                           log10_trans = FALSE,
                           negative = TRUE,
                           force_finite = FALSE,
                           includes_preds = FALSE,
                           suppress.messages = TRUE) {

  # Used variables such that there is no need to reference 'data'
  Chemical <- unique(data$Chemical)
  Species <- unique(data$Species)
  Time_trans <- data$Time_trans
  Dose <- data$Dose
  Route <- data$Route
  Media <- data$Media
  Conc <- data$Conc
  Conc_SD <- data$Conc_SD
  # Conc_trans <- data$Conc_trans # This is re-calculated here
  N_Subjects <- data$N_Subjects
  Detect <- data$Detect
  pLOQ <- data$pLOQ

  # Fix any NA in N_Subjects and Detect as we do in do_preprocess.pk()
  N_Subjects[is.na(N_Subjects)] <- 1.0
  Detect[is.na(Detect)] <- FALSE
  # Fix the modelfun_args and modelfun variables
  modelfun_args <- modelfun[["conc_fun_args"]]
  modelfun <- modelfun[["conc_fun"]]

  # combine parameters to be optimized and held constant
  params <- c(par, const_params)

  # Extract parameters whose names do not match 'sigma'
  # (that is, all the actual model parameters)
  sigma_index <- startsWith(names(params), prefix = "sigma")
  model.params <- params[!sigma_index]

  if (includes_preds) {
    stopifnot("data must have Conc_est column." = "Conc_est" %in% names(data))
    pred <- data$Conc_est
  } else {
    # assemble arguments for do.call
    if (!is.list(modelfun_args)) {
      modelfun_args <- as.list(modelfun_args)
    } else {
      # Need to evaluate some quoted arguments (CHEMICAL and SPECIES)
      # Need to evaluate it correctly in pk_model
      modelfun_args <- lapply(modelfun_args, rlang::eval_tidy, data = data)
    }
    these_args <- append(
      list(
        params = model.params,
        time = Time_trans,
        dose = Dose,
        route = Route,
        medium = Media
      ),
      modelfun_args
    )
    # get un-transformed predicted plasma concentration vs. time for the current parameter
    # values, by dose and route
    pred <- do.call(what = modelfun, args = these_args)
  }

  # Early Return
  # if model predicted any NA or Inf concentrations OR
  # any negative concentrations,
  # then parameters are infinitely unlikely
  if (!all(is.finite(pred)) || any(pred < 0)) {
    ll <- -Inf
    if (force_finite) ll <- -1 * sqrt(.Machine$double.xmax)
    if (negative) ll <- -1 * ll
    return(ll)
  }

  pred[which(pred < pLOQ)] <- pLOQ[which(pred < pLOQ)]

  # transform predictions and concentrations
  # dose normalization?
  if (dose_norm == TRUE) {
    pred_trans <- pred / Dose
    Conc_trans <- Conc / Dose
    Conc_SD_trans <- Conc_SD / Dose
  } else {
    pred_trans <- pred
    Conc_trans <- Conc
    Conc_SD_trans <- Conc_SD
  }

  # log10 transformation?
  if (log10_trans == TRUE) {
    pred_trans <- log10(pred_trans)
    Conc_trans <- log10(Conc_trans)
    Conc_SD_trans <- log10(Conc_SD_trans)
    ll_summary <- "dlnorm_summary"
    # maintain other transformations such as dose-scaling,
    # but undo log10 transformation
    conc_natural <- 10^(Conc_trans)
    conc_sd_natural <- 10^(Conc_SD_trans)
  } else {
    # if log10 transformation has *not* been applied
    ll_summary <- "dnorm_summary"
    conc_natural <- Conc_trans
    conc_sd_natural <- Conc_SD_trans
  }

  # residual error SDs
  # defined by data_sigma_group
  # Reversed conditional order for clarity
  if (!any(sigma_index)) {
    cli::cli_abort(c(
      "Could not find any parameters with 'sigma' in the name.",
      "i" = "Parameter name{?s} are: {names(params)}"
    ))
  }

  # There ARE sigma values
  # get sigma params
  sigma_params <- as.list(params[sigma_index])
  # match the study ID and assign each sigma to its corresponding study
  # if study ID doesn't match anything in the sigmas, then assign NA
  sigma_obs <- sigma_params[paste("sigma", data_sigma_group, sep = "_")]

  # sigma may be implicitly "known" if there is only one sigma value
  n_sigma_pars <- length(sigma_params)
  sigma_is_na <- vapply(sigma_obs, is.null, logical(1))
  sigma_obs[sigma_is_na] <- ifelse(n_sigma_pars == 1L, sigma_params, NA_real_)
  sigma_obs <- unlist(sigma_obs)

  # Refactoring to *more* vectorized version, there are six apparent conditions
  one_subject_bool = (N_Subjects == 1)
  # Put this in the context of either known or unknown sigmas
  deb_ks = one_subject_bool & Detect & !sigma_is_na # known sigma
  deb_us = one_subject_bool & Detect & sigma_is_na # unknown sigma
  ndb_ks = one_subject_bool & !Detect & !sigma_is_na
  ndb_us = one_subject_bool & Detect & sigma_is_na
  nsb_ks = !(one_subject_bool | sigma_is_na)
  nsb_us = !one_subject_bool & sigma_is_na

  ll_data_sigma <- sum(
    dnorm(
      x = Conc_trans[deb_ks],
      mean = pred_trans[deb_ks],
      sd = sigma_obs[deb_ks],
      log = TRUE
    ),
    pnorm(
      q = Conc_trans[ndb_ks],
      mean = pred_trans[ndb_ks],
      sd = sigma_obs[ndb_ks],
      log.p = TRUE
    ),
    do.call(
      ll_summary,
      list(
        mu = pred_trans[nsb_ks],
        sigma = sigma_obs[nsb_ks],
        x_mean = conc_natural[nsb_ks],
        x_sd = conc_sd_natural[nsb_ks],
        x_N = N_Subjects[nsb_ks],
        log = TRUE
      )
    )
  )

  # if data_sigma_group is NA -- then assume data are equally likely to come from
  # any of the existing distributions
  # data_sigma_group may be NA if we are calculating log-likelihood for new data,
  # i.e., data on which the model was not originally fitted.
  # compute log likelihoods for observations without sigmas
  if (any(sigma_is_na)) {
    if (isFALSE(suppress.messages)) {
      cli::cli_par()
      cli::cli_inform(
        c(
          paste("log_likelihood(): {sum(sigma_is_na)}, observation{?s}",
                "were not in any existing error-SD (sigma) group."),
          "i" = "They will be treated as equally likely to be in any of the existing error-SD groups."
        )
      )
      cli::cli_end()
    }

    ll_data_no_sigma <- vapply(
      unlist(sigma_params),
      function(this_sigma) {
        sum(
          dnorm(
            x = Conc_trans[deb_us],
            mean = pred_trans[deb_us],
            sd = sigma_obs[deb_us],
            log = TRUE
          ),
          pnorm(
            q = Conc_trans[ndb_us],
            mean = pred_trans[ndb_us],
            sd = sigma_obs[ndb_us],
            log.p = TRUE
          ),
          do.call(
            ll_summary,
            list(
              mu = pred_trans[nsb_us],
              sigma = sigma_obs[nsb_us],
              x_mean = conc_natural[nsb_us],
              x_sd = conc_sd_natural[nsb_us],
              x_N = N_Subjects[nsb_us],
              log = TRUE
            )
          )
        )
      },
      FUN.VALUE = numeric(1L)
    ) |> mean()
  } else {
    ll_data_no_sigma <- numeric(0L)
  }

  # sum log-likelihoods over observations
  ll <- sum(ll_data_sigma, ll_data_no_sigma)
  # do *not* remove NAs, because they mean this parameter combination is impossible!

  # If ll isn't finite,
  # just set it to -Inf to indicate that these parameters are infinitely unlikely
  if (!is.finite(ll)) ll <- -Inf

  # If user has selected to force return of a finite value,
  # e.g. as required by optimx with method 'L-BFGS-B',
  # then when log-likelihood is infinitely unlikely,
  # return a large negative number instead
  if (force_finite && !is.finite(ll)) {
    # now return sqrt of .Machine$double.xmax, not just .Machine$double.xmax
    # not taking sqrt seems to break L-BFGS-B sometimes
    ll <- -1 * sqrt(.Machine$double.xmax)
  }

  # to get negative log-likelihood (e.g. for minimization)
  if (negative) {
    ll <- -1 * ll
  }

  return(ll)
}
