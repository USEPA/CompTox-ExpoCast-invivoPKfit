#' Get upper bounds for estimating model parameters
#'
#' For a set of model parameters, get the upper bounds for the optimizer (if
#' parameter is to be estimated).
#'
#' The default `upper_default` `data.frame` is shown below in table format:
#'
#' | param_name     | upper_bound | upper_bound_msg |
#' | ---------------| ----------- | --------------- |
#' | kelim          | Inf        | Default         |
#' | Vdist          | Inf        | Default         |
#' | kgutabs        | Inf        | Default         |
#' | Fgutabs        | 1         | Default         |
#' | V1             | Inf           | Default         |
#' | k12            | Inf         | Default         |
#' | k21            | Inf         | Default         |
#' | Fgutabs_Vdist  | 1e4    | Default         |
#' | Fgutabs_V1     | 1e4    | Default         |
#' | sigma          | Inf         | Default         |
#'
#' #' Any parameters which will not be estimated from data (based on either the
#' variable `optimize_param` in `par_DF` if `par_DF` is provided, or the output
#' of [get_opt_params()] if `par_DF` is not provided) are assigned an upper bound
#' of `NA_real_`.
#'
#' @param fitdata A `data.frame`: the concentration-time-dose data to be used for
#'   fitting.
#' @param par_DF Optional: A data.frame as produced by [get_opt_params()], with a
#'   character variable `param_name` containing parameter names, and a logical
#'   variable `optimize_param` containing TRUE if parameter is to be fitted, and
#'   FALSE if parameter is to be held constant. Any other variables will be
#'   ignored. Default NULL, in which case it will be determined by calling
#'   [get_opt_params()].
#' @param model The name of the model whose parameters are to be estimated.
#'   Currently only "flat", "1compartment", or "2compartment" is supported.
#'   Ignored if `par_DF` is provided.
#' @param pool_sigma Logical: Whether to pool all data (estimate only one error
#'   standard deviation) or not (estimate separate error standard deviations for
#'   each reference). Default FALSE to estimate separate error SDs for each
#'   reference. (If `fitdata` only includes one reference, `pool_sigma` will
#'   have no effect.) Ignored if `par_DF` is provided.
#' @param upper_default A `data.frame` with three variables: `param_name`,
#'   giving the names of parameters; `upper_bound`, giving the default
#'   upper-bound values for each parameter; and `upper_bound_msg`, giving a
#'   message about the default upper bound values. See Details for default
#'   value.
#' @return A data.frame: `parDF` with additional variables `upper_bound`
#'   (numeric, containing the upper bound for each parameter) and
#'   `upper_bound_msg` (character, containing a brief message explaining how the
#'   upper-bound value was calculated).
#' @author Caroline Ring, John Wambaugh, Mitchell Teague
get_upper_bounds <- function(fitdata,
                             par_DF = NULL,
                             model,
                             pool_sigma = FALSE,
                             upper_default = data.frame(
                               param_name = c("kelim",
                                              "Vdist",
                                              "kgutabs",
                                              "Fgutabs",
                                              "V1",
                                              "k12",
                                              "k21",
                                              "Fgutabs_Vdist",
                                              "Fgutabs_V1",
                                              "Rblood2plasma",
                                              "sigma"),
                               upper_bound = c(1e6, #kelim
                                               1e6, #Vdist
                                               1e6, #kgutabs
                                               1, #Fgutabs
                                               1e6, #V1
                                               1e6, #k12
                                               1e6, #k21
                                               1e4, #Fgutabs_Vdist
                                               1e4, #Fgutabs_V1
                                               1e6, #Rblood2plasma
                                               1e6), #sigma
                               upper_bound_msg = "Default"
                             ),
                             sigma_from_data = TRUE,
                             fit_conc_dose = TRUE,
                             fit_log_conc = FALSE,
                             suppress.messages = FALSE){

  if(is.null(par_DF)){
  par_DF <- get_opt_params(model = model,
                           fitdata = fitdata,
                           pool_sigma = pool_sigma,
                           param_names = par_DF$param_name,
                           suppress.messages = suppress.messages)
  }
  rownames(par_DF) <- par_DF$param_name

  #Use the default upper bounds
  #this gets everything except sigma, which will not be in the model params
  par_DF <- merge(par_DF,
                  upper_default,
                  by = "param_name",
                  all.x = TRUE,
                  all.y = FALSE)


    #assign sigmas
    par_DF[grepl(x = par_DF$param_name,
                 pattern = "sigma"),
           "upper_bound"] <- upper_default[upper_default$param_name %in% "sigma",
                                           "upper_bound"]
    par_DF[grepl(x = par_DF$param_name,
                 pattern = "sigma"),
           "upper_bound_msg"] <- upper_default[upper_default$param_name %in% "sigma",
                                               "upper_bound_msg"]

    if(sigma_from_data %in% TRUE){
      #set an upper bound on sigma equal to the std dev of values in each study
      #on the grounds that it really should not be worse than that
      value_var <- ifelse(fit_conc_dose %in% TRUE, "Value_Dose", "Value")
      value_sd_var <- ifelse(fit_conc_dose %in% TRUE, "Value_SD_Dose", "Value_SD")
      sigma_names <- grep(x = par_DF$param_name,
                          pattern = "sigma",
                          value = TRUE)
      studies <- unique(fitdata$Study)
      if(length(sigma_names) > 1){ #if more than one study in this data set
#calculate each study SD, handling summary data properly
        data_sd <- sapply(studies,
                          function(x) {
                            foo <- as.list(fitdata[fitdata$Study %in% x,
                                           c(value_var,
                                             value_sd_var,
                                             "N_Subjects")])
                            names(foo) <- c("group_mean",
                                            "group_sd",
                                            "group_n")
                            do.call(combined_sd, args = c(foo,
                                                          list("log" = fit_log_conc)
                                                          )
                                    )
                            },
                          USE.NAMES = TRUE)
#name the study SDs after the sigmas
        names(data_sd) <- paste0("sigma_study_", studies)
        #assign each sigma
        for (sigma_name in names(data_sd)){
          if(is.finite(data_sd[sigma_name])){
            par_DF[par_DF$param_name %in% sigma_name,
                   "upper_bound"] <- 2*data_sd[sigma_name]
            par_DF[par_DF$param_name %in% sigma_name,
                   "upper_bound_msg"] <- paste("2*SD of", value_var,
                                               "for study",
                                               gsub(x = sigma_name,
                                                    pattern = "sigma_study_",
                                                    replacement = ""))
          }
        }
      }else{ #if only one sigma (one study, or all studies pooled)
        foo <- as.list(fitdata[,
                       c(value_var,
                         value_sd_var,
                         "N_Subjects")])
        names(foo) <- c("group_mean",
                        "group_sd",
                        "group_n")
        data_sd <- do.call(combined_sd,
                           args = c(foo,
                                    list("log" = fit_log_conc)
                           )
        )
        if(is.finite(data_sd)){
        par_DF[par_DF$param_name %in% "sigma", "upper_bound"] <- 2*data_sd
        par_DF[par_DF$param_name %in% "sigma",
               "upper_bound_msg"] <- paste("2*SD of", value_var,
                                           "for all studies in this data set")
        }
      }




    }

  #For anything not to be optimized, set its bounds to NA
  par_DF[!(par_DF$optimize_param %in% TRUE),
         c("upper_bound",
           "upper_bound_msg")] <- list(NA_real_,
                                       "optimize_param is not TRUE")

  #For anything not to be used, set its bounds to NA
  par_DF[!(par_DF$use_param %in% TRUE),
         c("upper_bound",
           "upper_bound_msg")] <- list(NA_real_,
                                       "use_param is not TRUE")

  return(par_DF)

}

