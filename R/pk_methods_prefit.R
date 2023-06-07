#' Pre-fitting
#'
#' Do pre-fit calculations and checks
#'
#' This function does the following:
#'
#' - Based on the error model in `stat_error_model` and the pre-processed data, determines the number of residual standard deviations ("sigmas") hyperparameters to be estimated.
#' - Determines which "sigma" hyperparameter corresponds to each observation in the data.
#' - Calculates lower/upper bounds and starting guesses for each "sigma" hyperparameter
#' - For each model in `stat_model`, calls its `params_fun`, the function that, based on the data, determines whether to optimize each model parameter, and calculates lower/upper bounds and starting guesses for each model parameter to be optimized. Only non-excluded observations are passed to each model's `params_fun`.
#'
#'
#' Lower bounds for each "sigma" hyperparameter are set to `sqrt(.Machine$double_eps)`.
#'
#' Upper bounds for each "sigma" hyperparameter are calculated as the standard
#' deviation of observations in the corresponding error SD group (see
#' [combined_sd()]). If the combined SD is non-finite or less than the sigma
#' lower bound, then the combined SD of all non-excluded data is substituted. If
#' that is still non-finite or less then the sigma lower bound, then a constant
#' value of 100 is substituted.
#'
#' The starting guess for each "sigma" hyperparameter is one-tenth of the upper bound.
#'
#' @param obj A `pk` object
#' @return The same `pk` object, but with a new element `prefit`, containing the
#'   results of pre-fit calculations and checks for each model and for the error
#'   model.
#' @export
#' @author Caroline Ring
prefit.pk <- function(obj){

  objname <- deparse(substitute(obj))
  status <- obj$status
  if(status >= status_prefit){
    warning(paste0(objname,
                   " current status is ",
                   status,
                   ". prefit() will reset its status to ",
                   status_prefit,
                   ". Any results from later workflow stages will be lost."))
  }

  #if preprocessing not already done, do it
  if(obj$status < status_preprocess){
    obj <- preprocess_data(obj)
  }

  if(obj$status < status_data_info){
    obj <- data_info(obj)
  }

  suppress.messages <- obj$settings_preprocess$suppress.messages

  data <- obj$data

  if(suppress.messages %in% FALSE){
    message(paste("prefit.pk(): Assigning error SD groups to all observations"))
  }

  #get the error model obj$stat_error_model, which defines the number of sigmas that will need to be optimized
  #count the number of unique combinations of vars in obj$stat_error_model$error_group
  unique_groups <- unique(rlang::eval_tidy(expr = obj$stat_error_model$error_group,
                                           data = data))

  #Assign a factor variable denoting sigma group to each observation in obj$data
  #This tells us which sigma applies to which observation
data_sigma_group <- interaction(
    lapply(
      obj$stat_error_model$error_group,
      function(x){
        rlang::eval_tidy(x, data = data)
      }
    )
  )

 #but set to NA for any excluded observations, and drop any levels that occur only for excluded observations
data_sigma_group[data$exclude %in% TRUE] <- NA_character_
data_sigma_group <- droplevels(data_sigma_group)

obj$prefit$stat_error_model$data_sigma_group <- data_sigma_group

if(suppress.messages %in% FALSE){
  message(paste("prefit.pk():",
                "Getting bounds and starting guesses for each error SD to be fitted"))
}
  #get bounds and starting points for each error sigma to be fitted
sigma_lower <- sqrt(.Machine$double.eps)
  sigma_DF <- data.frame(param_name = paste("sigma",
                                            levels(data_sigma_group),
                                            sep = "_"),
                         param_units = unique(data$Conc_trans.Units),
                         optimize_param = TRUE,
                         use_param = TRUE,
                         lower_bound = sigma_lower)
  rownames(sigma_DF) <- levels(data_sigma_group)

  #upper bounds: combined SD per group (except handle it if combined SD is zero)
  for(this_ds in rownames(sigma_DF)){
    DF_sub <- subset(data,
                     exclude %in% FALSE &
                       data_sigma_group %in% this_ds)
    if(nrow(DF_sub) > 0){
    sigma_upper <- combined_sd(
      group_mean = DF_sub$Conc_trans,
                                        group_sd = DF_sub$Conc_SD_trans,
                                        group_n = DF_sub$N_Subjects,
                                        unbiased = TRUE,
                                        na.rm = TRUE,
                                        log = FALSE)
    #if combined SD is non-finite or 0, then substitute with grand combined SD
    #of all non-excluded data
    if(!is.finite(sigma_upper) |
       sigma_upper <= sigma_lower){
      DF_sub <- subset(data,
                         exclude %in% FALSE)
      sigma_upper <- combined_sd(
        group_mean = DF_sub$Conc_trans,
        group_sd = DF_sub$Conc_SD_trans,
        group_n = DF_sub$N_Subjects,
        unbiased = TRUE,
        na.rm = TRUE,
        log = FALSE)
    }

    #if sigma_upper is still non-finite or 0, then impute 100
    if(!is.finite(sigma_upper) |
       sigma_upper <= sigma_lower){
      sigma_upper <- 100
    }
      sigma_DF[this_ds, "upper_bound"] <- sigma_upper
    }
  }



  #starting value = 0.1* upper bound
  sigma_DF$start <- 0.1 * sigma_DF$upper_bound

  #assign rownames to sigma_DF
  rownames(sigma_DF) <- sigma_DF$param_name

  #assign sigma_DF to the `pk` object
  obj$prefit$stat_error_model$sigma_DF <- sigma_DF

  n_sigma <- nrow(sigma_DF)
  #for each model to be fitted:
  for (this_model in names(obj$stat_model)){
    #get parameters to be optimized, bounds, and starting points
    #by evaluating params_fun for this stat_model
    #pass it only the non-excluded observations
    obj$prefit[[this_model]]$par_DF <- do.call(obj$stat_model[[this_model]]$params_fun,
                                                   args = c(list(subset(data, exclude %in% FALSE)),
                                                            obj$stat_model[[this_model]]$params_fun_args))
    #check whether there are enough observations to optimize the requested parameters plus sigmas
    #number of parameters to optimize
    n_par <- sum(obj$prefit[[this_model]]$par_DF$optimize_param)
    #number of detected, non-excluded observations
    n_detect <- sum(obj$data_info$data_summary_all$n_detect)
    obj$prefit[[this_model]]$fit_decision <- ifelse(
      n_detect <= (n_par + n_sigma),
      "abort",
      "continue"
    )

    obj$prefit[[this_model]]$fit_decision_reason <- ifelse(
      n_detect <= (n_par + n_sigma),
      paste0("Number of non-excluded detects (",
             n_detect,
             ") is less than or equal to number of parameters to optimize (",
             n_par, ") plus number of error SDs to optimize (",
             n_sigma,
             ")"),
      paste0("Number of non-excluded detects (",
             n_detect,
             ") is greater than number of parameters to optimize (",
             n_par, ") plus number of error SDs to optimize (",
             n_sigma,
             ")")
    )

    if(suppress.messages %in% FALSE){
    if(n_detect <= (n_par + n_sigma)){
      message(paste0("prefit.pk():",
                     "Model", this_model,
                     ": Fit will not be performed.",
                     "Number of non-excluded detects (",
                     n_detect,
                     ") is less than or equal to number of parameters to optimize (",
                     n_par, ") plus number of error SDs to optimize (",
                     n_sigma,
                     ")"))
    }
    }

  }

  obj$status <- status_prefit #prefit complete

  return(obj)

}
