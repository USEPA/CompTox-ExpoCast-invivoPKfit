#' Get starting values for fitting model parameters
#'
#' For a set of model parameters, get the starting guess for the optimizer (if
#' parameter is to be fitted) or the value at which to hold the parameter
#' constant (if parameter is not to be fitted).
#'
#' Model parameters are first assigned default starting values according to
#' `starts_default`.
#'
#' Then, the function checks to see if the package [httk] includes the necessary
#' TK data to parameterize a 1-compartment model for the chemical (determined by
#' `unique(fitdata$DTXSID)`. If so, then it assigns [httk]'s estimates for the
#' relevant model parameters.
#'
#' Finally, the function attempts to roughly estimate parameter starting values
#' from the data.
#'
#' The default `starts_default` `data.frame` is shown below in table form:
#'
#' | param_name     | start_value | start_value_msg |
#' | ---------------| ----------- | --------------- |
#' | A              | 1           | Default         |
#' | kelim          | 0.25        | Default         |
#' | Vdist          | 2.19        | Default         |
#' | kgutabs        | 0.25        | Default         |
#' | Fgutabs        | 0.5         | Default         |
#' | V1             | 1           | Default         |
#' | k12            | 0.2         | Default         |
#' | k21            | 0.5        | Default         |
#' | Fgutabs_Vdist  | 0.5/2.19    | Default         |
#' | Fgutabs_Vdist  | 0.5/2.19    | Default         |
#' | sigma          | 100         | Default         |
#'
#'
#' @param par_DF Optional: A data.frame as produced by [get_opt_params], with a
#'   character variable [param_name] containing parameter names, and a logical
#'   variable [optimize_param] containing TRUE if parameter is to be fitted, and
#'   FALSE if parameter is to be held constant. Any other variables will be
#'   ignored. Default NULL, in which case it will be determined by calling
#'   [get_opt_params].
#' @param model The name of the model whose parameters are to be estimated.
#' @param fitdata A data.frame: the concentration-time-dose data to be used for
#'   fitting.
#' @param starts_default A data.frame with variables `param_name` and
#'   `start_value`, defining all possible default parameter starting values
#'   across all models. See Details.
#' @param sigma_start_from_data Logical: TRUE means to override the starting
#'   value for "sigma" in `starts_default` and estimate a starting value as the
#'   median concentration in the data. FALSE means to keep the starting value
#'   for "sigma" in `starts_default`. Default TRUE.
#' @param suppress.messages Logical: Whether to suppress printing of messages.
#'
#' @return The data.frame `par_DF` with added variables `start_value`,
#'   containing each parameter's starting value, and `start_value_msg`,
#'   containing a brief message about how each starting value was calculated.
#' @author Caroline Ring, John Wambaugh, Mitchell Teague

get_starts <- function(par_DF = NULL,
                       model,
                       fitdata,
                       pool_sigma = FALSE,
                       starts_default = data.frame(
                         param_name = c("A",
                                        "kelim",
                                        "Vdist",
                                        "kgutabs",
                                        "Fgutabs",
                                        "V1",
                                        "k12",
                                        "k21",
                                        "Fgutabs_Vdist",
                                        "Fgutabs_V1",
                                        "sigma"),
                         start_value = c(1, #A
                                         0.25, #kelim
                                         5.56, #Vdist
                                         2.19, #kgutabs
                                         0.5, #Fgutabs
                                         1, #V1
                                         0.2, #k12
                                         0.5, #k21
                                         0.5/2.19, #Fgutabs_Vdist
                                         0.51/2.19, #Fgutabs_V1
                                         100), #sigma
                         start_value_msg = "Default"
                       ),
                       sigma_start_from_data = TRUE,
                       suppress.messages = FALSE
                       ){

  #if params not already specified, get them for this model
if(is.null(par_DF)){
    par_DF <- get_opt_params(model = model,
                             fitdata = fitdata,
                             pool_sigma = pool_sigma,
                             param_names = par_DF$param_name,
                             suppress.messages = suppress.messages)
}

  #get lower bounds if not already there
  if(!("lower_bound" %in% names(par_DF))){
    par_DF <- get_lower_bounds(par_DF = par_DF,
                               model = model,
                               pool_sigma = pool_sigma,
                               fitdata = fitdata,
                               suppress.messages = suppress.messages)
  }

  #get upper bounds if not already there
  if(!("upper_bound" %in% names(par_DF))){
    par_DF <- get_upper_bounds(par_DF = par_DF,
                               model = model,
                               pool_sigma = pool_sigma,
                               fitdata = fitdata,
                               suppress.messages = suppress.messages)
  }

  this.dtxsid <- unique(fitdata$DTXSID)
  if(length(this.dtxsid) > 1) stop("get_starts(): More than one DTXSID in data")

  this.species <- unique(fitdata$Species)
  if(length(this.species) > 1) stop("get_starts(): More than one species in data")

    #Use the default starting values to begin with
    par_DF <- merge(par_DF,
                    starts_default,
                    by = "param_name",
                    all.x = TRUE,
                    all.y = FALSE)

    rownames(par_DF) <- par_DF$param_name


    #if model is flat, take A to be the average concentration
  if(model %in% "flat"){
    par_DF["A", "start_value"] <- mean(fitdata$Value, na.rm = TRUE)
    par_DF["A", "start_value_msg"] <- "Mean concentration"
  }


    #if httk 1-comp model params exist, replace the defaults with these
    #and mark the source accordingly
    if(this.dtxsid %in% httk::get_cheminfo(info="dtxsid",
                                           model = "1compartment",
                                           suppress.messages = TRUE)){
      #use this_species if httk has data for it, otherwise use human
      httk_species <- ifelse(tolower(this.species) %in%
                               colnames(httk::physiology.data),
                             this.species,
                             "Human")
      #get httk parameters for this chemical & selected species
      suppressWarnings(httk_params <- httk::parameterize_1comp(dtxsid = this.dtxsid,
                                              default.to.human = TRUE,
                                              suppress.messages = TRUE,
                                              species = httk_species))
      if (model %in% "2compartment") {
        #set V1 to the httk-predicted Vdist
        #create a new "V1" parameter in httk_params for this purpose
        httk_params[["V1"]] <- httk_params[["Vdist"]]
      }
      #replace the parameters in par_DF whose names match those in httk_params
      #find matching names
      replace_names <- intersect(par_DF$param_name,
                                 names(httk_params))
      #do the replacement
      par_DF[match(replace_names,
                       par_DF$param_name,
                       nomatch = 0),
                 "start_value"] <- unlist(httk_params[replace_names])
      #and record the source of the new values (the httk function call)
      par_DF[match(replace_names,
                       par_DF$param_name,
                       nomatch = 0),
                 "start_value_msg"] <- paste0("httk::parameterize_1comp(",
                                                 "DTXSID = ", this.dtxsid, ", ",
                                                 "default.to.human = TRUE, ",
                                                 "suppress.messages = TRUE",
                                                 "species = ", httk_species,
                                                 ")")

    } #end if(this.dtxsid %in% httk::get_cheminfo(info="dtxsid"))

    #Try roughly estimating parameters from data

    #Fgutabs
      if(par_DF["Fgutabs", "optimize_param"] %in% TRUE &
         par_DF["Fgutabs", "use_param"] %in% TRUE){
    # Need both oral and iv data to get at Fgutabs
      if('po' %in% fitdata$Route &
         'iv' %in% fitdata$Route){
        iv.data <- subset(fitdata, Route == "iv" & !is.na(Value))
        iv.data <- subset(iv.data, Time == min(Time))
        iv.data <- subset(iv.data, Dose == max(Dose))
        oral.data <- subset(fitdata,Route == "po" & !is.na(Value))
        oral.data <- subset(oral.data, Value == max(Value))[1, ]
        Fgutabs_tmp <- oral.data$Value / mean(iv.data$Value) *
          mean(iv.data$Dose) / mean(oral.data$Dose)
        #if the estimated value is valid, then use it and record the source
        if (is.finite(Fgutabs_tmp) &
            Fgutabs_tmp > par_DF["Fgutabs", "lower_bound"] &
            Fgutabs_tmp < par_DF["Fgutabs", "upper_bound"]){
          par_DF["Fgutabs",
                     "start_value"] <- Fgutabs_tmp
          par_DF["Fgutabs",
                     "start_value_msg"] <- "Estimated from oral and IV data"
        }else{
          par_DF["Fgutabs",
                     "start_value"] <- 0.5
          par_DF["Fgutabs",
                     "start_value_msg"] <- paste0("Fgutabs assumed = 0.5, ",
                     "because estimation from oral and IV data ",
                     "produced either NaN or a value outside [0,1]")
        }
      }else{
        warn_msg <- "get_starts(): Cannot optimize Fgutabs without both oral and IV data"
        warning(warn_msg)
        par_DF["Fgutabs", "start_value"] <- NA_real_
        par_DF["Fgutabs", "start_value_msg"] <- warn_msg
      } #end if('po' %in% fitdata$Route & 'iv' %in% fitdata$Route) / else
      } #end if(par_DF["Fgutabs", "optimize_parambo"] %in% TRUE)


      #kelim
      if(par_DF["kelim", "optimize_param"] %in% TRUE){
        #make elim data
        elim_data <- make_elim_data(fitdata = fitdata)

        #take only the later times
        time_split <- median(elim_data$Time)
        elim_late <- subset(elim_data,
                            Time >= time_split)

    #linear regression of log conc/dose vs. time
    elim_lm <- stats::lm(log(ValueDose) ~ Time,
                         data=elim_data)

    elim_coeff <- coef(elim_lm)

    #take negative slope
    kelim_tmp <- -elim_coeff[2]

    #if kelim_tmp is not NA, then use it.
    if(is.finite(kelim_tmp)){
      if(kelim_tmp > par_DF["kelim", "lower_bound"]){
      par_DF[par_DF$param_name %in% "kelim",
                 "start_value"] <- kelim_tmp
      par_DF[par_DF$param_name %in% "kelim",
                 "start_value_msg"] <- paste0("Estimated from linear regression on ",
                                                 unique(elim_data$Route),
                                                 " data")
      }else{
        par_DF[par_DF$param_name %in% "kelim",
                   "start_value"] <- par_DF["kelim", "lower_bound"]
        par_DF[par_DF$param_name %in% "kelim",
                   "start_value_msg"] <- paste0("Linear regression on ",
                                                unique(elim_data$Route),
                                                   " data produced NA or < ",
                                                "lower bound",
                                                   "; setting to lower bound")
      }
    }
      }


      #### vdist or V1
  if(any(par_DF[c("Vdist","V1"),
              "optimize_param"] %in% TRUE)){
    #make elim data
    elim_data <- make_elim_data(fitdata = fitdata)
    #take the avg conc/dose at the earliest timepoint
    Vd_data <- subset(elim_data, Time == min(Time))
    Vd_tmp <- 1 / mean(Vd_data$ValueDose)
    #select which parameter to assign: Vdist or V1
    Vd_name <- intersect(par_DF$param_name,
                         c("Vdist", "V1"))
    if(Vd_tmp <= par_DF[par_DF$param_name %in% Vd_name,
                        "upper_bound"] &
       (Vd_tmp >= par_DF[par_DF$param_name %in% Vd_name,
                         "lower_bound"])){
      #assign starting value
      par_DF[Vd_name,
             "start_value"] <- Vd_tmp
      #record source
      par_DF[Vd_name,
             "start_value_msg"] <- paste0("Estimated from ",
                                          unique(Vd_data$Route),
                                          " data")
    }
  }

      #otherwise, if we are estimating Fgutabs/V, use the Vd data that way
      if(any(par_DF[c("Fgutabs_Vdist","Fgutabs_V1"),
                    "optimize_param"] %in% TRUE)){
        #make elim data
        elim_data <- make_elim_data(fitdata = fitdata)
        #use the same data as for kelim (transformed in the same way)
        #take the avg conc/dose at the earliest timepoint
        FV_data <- subset(elim_data, Time == min(Time))
        FV_tmp <- mean(FV_data$ValueDose)
        #select which parameter to assign: Vdist or V1
        FV_name <- intersect(par_DF$param_name,
                             c("Fgutabs_Vdist", "Fgutabs_V1"))
        if(FV_tmp <= par_DF[par_DF$param_name %in% FV_name,
                            "upper_bound"] &
           (FV_tmp >= par_DF[par_DF$param_name %in% FV_name,
                             "lower_bound"])){
          #assign starting value
          par_DF[FV_name,
                 "start_value"] <- FV_tmp
          #record source
          par_DF[FV_name,
                 "start_value_msg"] <- paste0("Estimated from ",
                                              unique(FV_data$Route),
                                              " data")
        }
      }


      #k21, kel, k12
      if(all(par_DF[c("kelim",
                  "k12",
                  "k21"),
                  "optimize_param"] %in% TRUE)){
        #make elim data
        elim_data <- make_elim_data(fitdata = fitdata)

        if (length(unique(elim_data$Time)) > 4) {
          #if we have enough time points to separate early and late:
          #use the elim_data prepared earlier for kelim
          #get the method of residuals:
          #first take late-phase data: the last half of time
          time_split <- median(elim_data$Time)
          elim_late <- subset(elim_data,
                              Time >= time_split)
          beta_lm <- stats::lm(log(Value) ~ Time,
                                data =  elim_late
          )

          beta <- -coef(beta_lm)[2]
          B <- coef(beta_lm)[1]

          #then take early-phase data
          elim_early <- subset(elim_data,
                               Time <= time_split)
          #subtract conc - prediction of beta_lm
          elim_early$Value_beta <- predict(beta_lm,
                                           newdata = data.frame(
                                             Time = elim_early$Time
                                             )
                                           )
          #residuals: if not positive, make them small, because we have to log-transform
          elim_early$Value_resid <- pmax(elim_early$Value -
                                           elim_early$Value_beta,
                                         1e-8)

          #regression of residuals
          alpha_lm <- stats::lm(log(Value_resid) ~ Time,
                                data = elim_early)

          alpha <- -coef(alpha_lm)[2]
          A <- coef(alpha_lm)[1]

          ABsum <- A+B

          tmp_msg <- paste("A, B, alpha, beta estimated using method of residuals.",
                           "Early and late times split at median time value in data.")

        }else{
          alpha <- 2
          beta <- 0.5
          Vd_data <- subset(elim_data, Time == min(Time))
          ABsum <-  mean(Vd_data$ValueDose)
          A <- 0.75 * ABsum
          B <- 0.25 * ABsum

          tmp_msg <- paste("Could not estimate with only ",
                                length(unique(elim_data$Time)),
                                "time points. Assuming",
                           "alpha = 2, A = 0.75*mean(Value/Dose) at min time,",
                           "beta = 0.5, B = 0.25*mean(Value/Dose) at min time.")
        }
        k21 <-( A*beta + B*alpha)/ABsum
        kelim <- (alpha*beta)/k21
        k12 <- alpha + beta - k21 - kel

        if(k21 >= par_DF["k21", "lower_bound"] &
           k21 <= par_DF["k21", "upper_bound"]){
          par_DF[ "k21",
                 c("start_value",
                   "start_value_msg")] <- list(k21, tmp_msg)
        } #else retain default value

        if(kelim >= par_DF["kelim", "lower_bound"] &
           kelim <= par_DF["kelim", "upper_bound"]){
          par_DF[ "kelim",
                  c("start_value",
                    "start_value_msg")] <- list(kelim, tmp_msg)
        }#else retain default value

        if(k12 >= par_DF["k12", "lower_bound"] &
           k12 <= par_DF["k12", "upper_bound"]){
          par_DF[ "k12",
                  c("start_value",
                    "start_value_msg")] <- list(k12, tmp_msg)
        }#else retain default value
      }

      if(sigma_start_from_data %in% TRUE){
      # Initialize the standard deviations to the median concentration
  TRYSIGMA <- stats::median(fitdata$Value, na.rm = TRUE)
  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         c("start_value",
           "start_value_msg")] <- list(TRYSIGMA,
                                       paste0("median conc"))
      }else{
        par_DF[grepl(x = par_DF$param_name,
                     pattern = "sigma"),
               c("start_value",
                 "start_value_msg")] <- starts_default[
                   starts_default$param_name %in% "sigma",
                   c("start_value",
                     "start_value_msg")]
      }

  #For anything with use_param == FALSE, set the starting value to NA
  par_DF[!(par_DF$use_param %in% TRUE),
         c("start_value",
           "start_value_msg")] <- list(NA_real_,
                                    "use_param is not TRUE")


  return(par_DF)
}
