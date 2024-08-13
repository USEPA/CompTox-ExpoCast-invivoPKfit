#' Do pre-fitting
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
#' @param ... Additional arguments. Not in use.
#' @return The same `pk` object, but with a new element `prefit`, containing the
#'   results of pre-fit calculations and checks for each model and for the error
#'   model.
#' @export
#' @author Caroline Ring
do_prefit.pk <- function(obj,
                         ...){

  objname <- deparse(substitute(obj))
  status <- obj$status
  if(status >= status_prefit){
    warning(paste0(objname,
                   " current status is ",
                   status,
                   ". do_prefit.pk() will reset its status to ",
                   status_prefit,
                   ". Any results from later workflow stages will be lost."))
  }

  #if preprocessing not already done, do it
  if(obj$status < status_preprocess){
    obj <- do_preprocess(obj)
  }

  if(obj$status < status_data_info){
    obj <- do_data_info(obj)
  }

  suppress.messages <- obj$settings_preprocess$suppress.messages

  data <- obj$data

  if(suppress.messages %in% FALSE){
    message(paste("do_prefit.pk(): Assigning error SD groups to all observations"))
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

#add data_sigma_group to data
obj$data$data_sigma_group <- data_sigma_group
data <- get_data(obj)

#get bounds and starting points for each error sigma to be fitted

if(suppress.messages %in% FALSE){
  message(paste("do_prefit.pk():",
                "Getting bounds and starting guesses for each error SD to be fitted"))
}

# Set a value to square root of lowest possible value 'x' where 1+x != 1
sigma_lower <- sqrt(.Machine$double.eps)

# Add the data_sigma_group and filter out all excluded values
sigma_DF <- data %>%
  dplyr::mutate(data_sigma_group = data_sigma_group) %>%
  dplyr::filter(exclude %in% FALSE) %>%
  #temporarily undo log10-trans, if it has been used
  dplyr::mutate(Conc_tmp = dplyr::if_else(rep(obj$scales$conc$log10_trans %in% TRUE,
                                              NROW(Conc_trans)),
                                  10^Conc_trans,
                                  Conc_trans),
                Conc_SD_tmp = dplyr::if_else(rep(obj$scales$conc$log10_trans %in% TRUE,
                                                 NROW(Conc_trans)),
                                     10^Conc_SD_trans,
                                     Conc_SD_trans),
                Conc_tmp.Units = dplyr::if_else(rep(obj$scales$conc$log10_trans %in% TRUE,
                                                    NROW(Conc_trans)),
                                                gsub(pattern = "\\)$",
                                                     replacement = "",
                                                     x = gsub(pattern = "^log10\\(",
                                                     replacement = "",
                                                     x = Conc_trans.Units)),
                                                Conc_trans.Units))

# Set values for sigma upper/lower-bounds and start
sigma_DF <- do.call(dplyr::group_by,
                    args = c(list(sigma_DF),
                             obj$stat_error_model$error_group)) %>%
  dplyr::summarise(param_name = paste("sigma",
                                      unique(data_sigma_group),
                                      sep = "_"),
                   param_units = unique(Conc_tmp.Units),
                   optimize_param = TRUE,
                   use_param = TRUE,
                   lower_bound = sigma_lower,
                   upper_bound = combined_sd(
                     group_mean = Conc_tmp,
                     group_sd = Conc_SD_tmp,
                     group_n = N_Subjects,
                     unbiased = TRUE,
                     na.rm = TRUE,
                     log10 = obj$scales$conc$log10_trans),
                   start = 0.1 * upper_bound) %>%
  as.data.frame()

  #assign rownames to sigma_DF
  rownames(sigma_DF) <- sigma_DF$param_name

  #assign sigma_DF to the `pk` object
  obj$prefit$stat_error_model$sigma_DF <- sigma_DF


  if(suppress.messages %in% FALSE){
    message(paste("do_prefit.pk():",
                  "Getting bounds and starting guesses for all model parameters to be fitted"))
  }

  #for each model to be fitted:
  par_DF_out <- sapply(names(obj$stat_model),
         function(this_model){
    #get parameters to be optimized, bounds, and starting points
    #by evaluating params_fun for this stat_model
    #pass it only the non-excluded observations
    par_DF <- data %>%
      dplyr::filter(exclude %in% FALSE)

    par_DF <- dplyr::group_by(par_DF,
                              !!!obj$data_group) %>%
      dplyr::reframe(do.call(obj$stat_model[[this_model]]$params_fun,
                               args = c(list(dplyr::cur_data_all()),
                                        obj$stat_model[[this_model]]$params_fun_args)
                               )
      ) %>% as.data.frame()

     par_DF
         },
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  par_DF_out <- do.call(dplyr::bind_rows,
                        c(par_DF_out,
                        list(.id = "model")))

  fit_check_out <-  sapply(names(obj$stat_model),
                           function(this_model){
    #check whether there are enough observations to optimize the requested parameters plus sigmas
    #number of parameters to optimize
                             if(suppress.messages %in% FALSE){
                               message(paste("do_prefit.pk():",
                                             "Checking whether sufficient observations to fit models"))
                             }

                             n_par_DF <- par_DF_out %>%
                               dplyr::filter(model %in% this_model) %>%
                               dplyr::group_by(!!!obj$data_group) %>%
                               dplyr::summarise(n_par = sum(optimize_param))

                             n_sigma_DF <- sigma_DF %>%
                               dplyr::group_by(!!!obj$data_group) %>%
                               dplyr::summarise(n_sigma = sum(optimize_param))


                             n_detect_DF <- get_data_summary(obj) %>%
                               dplyr::group_by(!!!obj$data_group) %>%
                               dplyr::summarise(n_detect = sum(n_detect))


                             #merge all of these together
                             fit_check_DF <- dplyr::inner_join(
                               dplyr::inner_join(n_par_DF,
                                                 n_sigma_DF,
                                                 by = sapply(obj$data_group,
                                                             rlang::as_label)),
                               n_detect_DF,
                               by = sapply(obj$data_group,
                                           rlang::as_label)
                             )

                             #get fit decision & reasoning
                             fit_check_DF <- fit_check_DF %>%
                               dplyr::mutate(n_par_opt = n_par + n_sigma,
                                             fit_decision = ifelse(n_par_opt < n_detect,
                                                                   "continue",
                                                                   "abort"),
                                             fit_reason = ifelse(n_par_opt < n_detect,
                                                                 "Number of parameters to estimate is less than number of non-excluded detected observations",
                                                                 "Number of parameters to estimate is greater than or equal to number of non-excluded detected observations")) %>%
                               as.data.frame()

                             fit_check_DF
                           },
    simplify = FALSE,
    USE.NAMES = TRUE)

  fit_check_out <- do.call(dplyr::bind_rows,
                           c(fit_check_out,
                             list(.id = "model")))

  obj$prefit$par_DF <- par_DF_out
  obj$prefit$fit_check <- fit_check_out

  obj$status <- status_prefit #prefit complete

  return(obj)

}
