#' Do pre-fit calculations and checks
#'
#' @param obj A `pk` object
#' @return The same `pk` object, but with additional elements added to each item
#'   in `models`, containing the results of pre-fit calculations and checks for
#'   each model.
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

  #get the error model obj$stat_error_model, which defines the number of sigmas that will need to be optimized
  #count the number of unique combinations of vars in obj$stat_error_model$error_group
  unique_groups <- unique(rlang::eval_tidy(expr = obj$stat_error_model$error_group,
                                           data = obj$data))
  n_sigma <- nrow(unique_groups)

  obj$prefit$stat_error_model$n_sigma <- n_sigma

  #Assign a factor variable denoting sigma group to each observation in obj$data
  #This tells us which sigma applies to which observation
 obj$prefit$stat_error_model$data_sigma_group <- interaction(
    lapply(
      obj$stat_error_model$error_group,
      function(x){
        rlang::eval_tidy(x, data = obj$data)
      }
    )
  )

  #get bounds and starting points for each error sigma to be fitted
  sigma_DF <- data.frame(param_name = paste("sigma",
                                            levels(obj$prefit$stat_error_model$data_sigma_group),
                                            sep = "_"),
                         param_units = unique(obj$data$Conc_trans.Units),
                         optimize_param = TRUE,
                         use_param = TRUE,
                         lower_bound = .Machine$double.eps)

  #get upper bound: standard deviation of the transformed concentration
  #(Conc_trans) in each group
  sigma_DF$upper_bound <- tapply(X = obj$data$Conc_trans,
                                 INDEX =obj$prefit$stat_error_model$data_sigma_group,
                                 FUN = sd,
                                 na.rm = TRUE,
                                 simplify = TRUE)

  #get starting value for sigma: say, 0.5 of the upper bound
  sigma_DF$start <- 0.5*sigma_DF$upper_bound

  #assign rownames to sigma_DF
  rownames(sigma_DF) <- sigma_DF$param_name

  #assign sigma_DF to the `pk` object
  obj$prefit$stat_error_model$sigma_DF <- sigma_DF

  n_sigma <- nrow(sigma_DF)
  #for each model to be fitted:
  for (this_model in names(obj$stat_model)){
    #get parameters to be optimized, bounds, and starting points
    #by evaluating params_fun for this stat_model
    obj$prefit[[this_model]]$par_DF <- do.call(obj$stat_model[[this_model]]$params_fun,
                                                   args = c(list(obj$data),
                                                            obj$stat_model[[this_model]]$params_fun_args))
    #check whether there are enough observations to optimize the requested parameters plus sigmas
    #number of parameters to optimize
    n_par <- sum(obj$prefit[[this_model]]$par_DF$optimize_param)
    #number of detected observations
    n_detect <- sum(obj$data$Detect)
    obj$prefit[[this_model]]$fit_decision <- ifelse(
      n_detect <= (n_par + n_sigma),
      "abort",
      "continue"
    )

    obj$prefit[[this_model]]$fit_decision_reason <- ifelse(
      n_detect <= (n_par + n_sigma),
      paste0("Number of detects (",
             n_detect,
             ") is less than or equal to number of parameters to optimize (",
             n_par, ") plus number of error SDs to optimize (",
             n_sigma,
             ")"),
      paste0("Number of detects (",
             n_detect,
             ") is greater than number of parameters to optimize (",
             n_par, ") plus number of error SDs to optimize (",
             n_sigma,
             ")")
    )

  }

  obj$status <- status_prefit #prefit complete

  return(obj)

}
