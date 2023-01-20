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
#' @param Fgutabs_Vdist_from_species Logical: TRUE to estimate the upper bound
#'   of `Fgutabs_Vdist` or `Fgutabs_V1` using the upper bound of `Fgutabs`
#'   specified in `upper_default`, divided by a species-specific lower bound on
#'   the volume of distribution, based on species-specific physiological plasma
#'   volumes. FALSE to use the upper bound for `Fgutabs_Vdist` or `Fgutabs_V1`
#'   specified in `upper_default`.
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
                               upper_bound = c(Inf, #kelim
                                               Inf, #Vdist
                                               Inf, #kgutabs
                                               1, #Fgutabs
                                               Inf, #V1
                                               Inf, #k12
                                               Inf, #k21
                                               1e4, #Fgutabs_Vdist
                                               1e4, #Fgutabs_V1
                                               Inf, #Rblood2plasma
                                               Inf), #sigma
                               upper_bound_msg = "Default"
                             ),
                             Fgutabs_Vdist_from_species = FALSE,
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

    if(Fgutabs_Vdist_from_species %in% TRUE &
       model %in% c("1compartment", "2compartment")){
      #For Vdist or V1: Set the theoretical lower bound to something on the order of
      #the species-specific total plasma volume, pulled from httk physiology.data
      phys <- httk::physiology.data
      species <- unique(tolower(fitdata$Species)) #there should be only one species
      #if httk::physiology.data has data on this species:
      if(any(tolower(names(phys)) %in%  species)){
        #get the column number corresponding to this species
        species_col <- which(tolower(names(phys)) %in%
                               species)
        #pull the average plasma volume in mL/kg
        plasma_vol_mL_kg <- phys[phys$Parameter %in% "Plasma Volume",
                                 species_col]
        #pull the average body weight in kg
        body_wt_kg <- phys[phys$Parameter %in% "Average BW",
                           species_col]
        #convert mL/kg to L plasma
        plasma_vol_L <- plasma_vol_mL_kg * body_wt_kg * 1e-3
        Vdist_lower <- plasma_vol_L/10
        Fgutabs_upper <- par_DF[par_DF$param_name %in% "Fgutabs",
                                "upper_bound"]
        #assuming this gives us a valid answer:
        if(is.finite(Vdist_lower)){
          #theoretical upper bound on Fgutabs/Vdist is 1/Vdist_lower
          par_DF[par_DF$param_name %in%
                   c("Fgutabs_Vdist", "Fgutabs_V1"),
                 "upper_bound"] <- Fgutabs_upper/Vdist_lower
          #add message explaining lower bound & source
          par_DF[par_DF$param_name %in%
                   c("Fgutabs_Vdist", "Fgutabs_V1"),
                 "upper_bound_msg"] <- paste0("[Upper bound of Fgutabs,",
                                              signif(Fgutabs_upper, digits = 3),
                 "]/ [lower bound of Vdist, 0.1 * ",
                                              "average total plasma volume for species (",
                                              signif(plasma_vol_L, digits = 3),
                                              " L), from httk::physiology.data]")
        } #end if(is.finite(plasma_vol_L))
      } #end if(any(tolower(names(phys)) %in%  species))
    } #end if if(Fgutabs_Vdist_from_species %in% TRUE)

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
                                            "group_N")
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
                   "upper_bound"] <- data_sd[sigma_name]
            par_DF[par_DF$param_name %in% sigma_name,
                   "upper_bound_msg"] <- paste("SD of", value_var,
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
                        "group_N")
        data_sd <- do.call(combined_sd,
                           args = c(foo,
                                    list("log" = fit_log_conc)
                           )
        )
        if(is.finite(data_sd)){
        par_DF[par_DF$param_name %in% "sigma", "upper_bound"] <- data_sd
        par_DF[par_DF$param_name %in% "sigma",
               "upper_bound_msg"] <- paste("SD of", value_var,
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

#' Combined standard deviation
#'
#' Given mean, standard deviation, and N for some set of groups, calculate the
#' combined standard deviation. Note that the groups may not overlap.
#'
#' @param group_mean A numeric vector of group means.
#' @param group_sd A numeric vector of group SDs.
#' @param group_N A numeric vector of group sizes.
#' @param unbiased Logical. If TRUE, then `group_sd` is assumed to be the
#'   unbiased estimator of population standard deviation (i.e. using `n-1` in
#'   the denominator -- the way that `stats::sd()` calculates it), and the
#'   returned combined SD is also the unbiased estimator of the combined
#'   population SD. If FALSE, then `group_sd` is assumed to be the biased
#'   estimator (using `n` in the denominator), and the returned value is also
#'   the biased estimator of the combined population SD.
#' @param na.rm Logical. If TRUE (default), then any groups where mean, SD, *or* N
#' were NA will be dropped. If FALSE, they will be retained (and the result will
#' be NA).
#' @return Numeric: the standard deviation of the combined population (i.e. if
#'   all the groups were concatenated into one large group).
#' @author Caroline Ring
combined_sd <- function(group_mean,
                       group_sd,
                       group_N,
                       unbiased = TRUE,
                       na.rm = TRUE,
                       log = FALSE){

  x_len <- c("group_mean" = length(group_mean),
             "group_sd" = length(group_sd),
             "group_N" = length(group_N))

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
  which_na <- is.na(group_mean) | is.na(group_sd) | is.na(group_N)
  group_mean <- group_mean[!which_na]
  group_sd <- group_sd[!which_na]
  group_N <- group_N[!which_na]
  }

  grand_mean <- sum(group_N*group_mean)/sum(group_N)

  #if all N = 1, then just take regular standard deviation
  if(all(group_N %in% 1)){
    grand_sd <- sd(group_mean)

    grand_N <- sum(group_N)
    if(unbiased %in% FALSE){
      #convert unbiased SD to biased SD
      grand_sd <- grand_sd * sqrt((grand_N-1)/grand_N)
    }
  }else{ #if not all N = 1
  if(unbiased %in% TRUE){
    #convert unbiased group SDs to biased group SDs
    group_sd <- group_sd *
      sqrt((group_N[which_na]-1)/group_N)
  }

  grand_var <- (sum(group_N*group_sd^2) +
      sum(group_N*(group_mean - grand_mean)^2))/
    (sum(group_N))

  grand_sd <- sqrt(grand_var) #biased


  if(unbiased %in% TRUE){
    grand_N <- sum(group_N)
    #convert biased grand SD to unbiased grand SD
grand_sd <- grand_sd * sqrt(grand_N/(grand_N - 1))
  }
  }

  if(log %in% TRUE){
    #convert to log-scale combined SD
    grand_sd <- sqrt(log(1 + grand_sd^2 / grand_mean^2))
  }

  return(grand_sd)
}
