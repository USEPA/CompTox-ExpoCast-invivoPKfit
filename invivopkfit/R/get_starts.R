#' Get starting values for fitting model parameters
#'
#' For a set of model parameters, get the starting guess for the optimizer (if
#' parameter is to be fitted) or the value at which to hold the parameter
#' constant (if parameter is not to be fitted).
#'
#' Starting values are estimated in the following way.
#'
#' # Default starting values
#'
#' Model parameters are first assigned default starting values according to
#' `starts_default`.
#'
#'
#' The default `starts_default` `data.frame` is shown below in table format:
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
#' | k21            | 0.5         | Default         |
#' | Fgutabs_Vdist  | 0.5/2.19    | Default         |
#' | Fgutabs_V1     | 0.5/2.19    | Default         |
#' | sigma          | 1         | Default         |
#'
#'
#' # `httk` 1-compartment model parameterization starting values
#'
#' Then, for any parameters named in `start_from_httk`, the function checks,
#' using [httk::get_cheminfo()], to see if the package [httk] includes the
#' necessary TK data to parameterize a 1-compartment model for the chemical
#' (determined by `unique(fitdata$DTXSID)`). If so, then it uses
#' [httk::parameterize_1comp()] to get [httk]'s parameter estimates, and uses
#' them as the starting values for any parameters named in `start_from_httk`.
#' Parameters that can be estimated using [httk::parameterize_1comp()] include
#' `Vdist`, `kelim`, `kgutabs`, and `Fgutabs`. If `model == '2compartment'`,
#' then the `httk` estimate of `Vdist` will be used as the starting value of
#' `V1`.
#'
#' # Data-derived starting values
#'
#' Finally, for any parameters named in `start_from_data`, the function attempts
#' to roughly estimate parameter starting values from the data, with methods
#' depending on the model selected in `model` and the dosing routes contained in
#' the data in `fitdata`.
#'
#' ## Flat model
#'
#' A starting value for `A` is estimated as the median concentration in `fitdata`.
#'
#' ## 1-compartment model
#'
#' `fitdata` is divided into two parts: IV dosing data and PO (oral) dosing
#' data. (If only one of the two routes is available, then only the available
#' data is used.) Control data are removed (all observations where `Dose == 0`).
#' Concentrations are normalized by Dose, then log-transformed.
#'
#' ### IV dosing data available
#'
#' `Vdist` and `kelim` are estimated using linear regression on
#' `log(Value/Dose)` vs. `Time` from the IV data. `kelim` is the negative slope
#' of this linear regression; `Vdist` is the intercept, exponentiated to bring
#' it back to the natural scale.
#'
#' ### PO dosing data available
#'
#' `kelim`, `kgutabs`, and `Fgutabs_Vidst` (the ratio of `Fgutabs` to `Vdist`)
#' are estimated using the method of residuals on `log(Value/Dose)` vs. `Time`
#' from the PO data.
#'
#' Specifically, the PO data are divided into two parts: "early" and "late"
#' phases, using the time of peak concentration as the dividing line.
#'
#' First, a linear regression is performed on the "late" data. The negative
#' slope of the "late" linear regression gives an estimate for `kelim`. The
#' intercept of the "late" regression is exponentiated and recorded as `A`.
#'
#' Then, the "early" data are predicted using the "late" linear
#' regression, and the residuals are calculated. A second linear regression is
#' performed on the "early" residuals. The negative slope of this second linear
#' regression provides an estimate for `kgutabs`. Then, the intercept of the
#' late regression, `A`, can be used with the estimates for `kelim` and
#' `kgutabs` to calculate an estimate for `Fgutabs_Vdist`:
#'
#' -  `Fgutabs_Vdist <- A * (kgutab - kelim) / kgutabs`
#'
#' #### Both PO and IV dosing data available
#'
#' The final starting value for `kelim` is taken as the average of the IV-based
#' estimate and the PO-based estimate (with any NAs removed).
#'
#' The estimates for `Vdist` (from the IV data) and `Fgutabs_Vdist` (from the PO
#' data) are multiplied to produce an estimate for `Fgutabs`.
#'
#' ## 2-compartment model
#'
#' `fitdata` is divided into two parts: IV dosing data and PO (oral) dosing
#' data. (If only one of the two routes is available, then only the available
#' data is used.) Control data are removed (all observations where `Dose == 0`).
#' Concentrations are normalized by Dose, then log-transformed.
#'
#' ### IV dosing data available
#'
#' `V1`, `kelim`, `k12` and `k21` are estimated using the method of residuals on
#' `log(Value/Dose)` vs. `Time` from the PO data. Specifically, the data are
#' divided into "early" and "late" parts, with the dividing line being the
#' "elbow point" of the data (using [akmedoids::elbow_point()]).
#'
#' First, a
#' linear regression is performed on the "late" data. The negative slope of the
#' "late" linear regression is recorded as `beta`, and the exponentiated
#' intercept as `B`.
#'
#' Then the "early" data are predicted using the "late" linear
#' regression, and the residuals are calculated. A second linear regression is
#' performed on the "early" residuals. The negative slope of this second linear
#' regression is recorded as `alpha`, and the exponentiated intercept as `A`.
#' Then, the parameter starting values are calculated as follows (see
#' https://www.boomer.org/c/p4/c19/c1902.php):
#'
#' -`k21 <- (A * beta + B*alpha)/(A + B)`
#' -`kelim <- (alpha * beta)/k21`
#' -`k12 <- alpha + beta - k21 - kel`
#' -`V1 <- (alpha - k21)/(A * (alpha - beta))`
#'
#' In case there are not enough data available in each phase to perform regression,
#'
#' ### PO dosing data available
#'
#' `kelim`, `k12`, `k21` `kgutabs`, and `Fgutabs_V1` (the ratio of `Fgutabs` to
#' `V1`) are estimated using the method of residuals on `log(Value/Dose)` vs.
#' `Time` from the PO data. Specifically, the data are divided into three parts:
#' "absorption" phase, "early" phase, and "late" phase. The time of peak
#' concentration is taken as the end of the "absorption" phase data. Then, for
#' all data after the time of peak concentration, the "elbow" point is found
#' (using [akmedoids::elbow_point()]) and used as the dividing line between
#' early and late-phase data.
#'
#' First, a linear regression is performed on the "late" data. The negative
#' slope of the "late" linear regression is recorded as `beta`, and the
#' exponentiated intercept as `B`.
#'
#' Then the "early"-phase data are predicted using the "late" linear regression,
#' and the residuals are calculated. A second linear regression is performed on
#' the "early" residuals. The negative slope of this "early"-residuals linear
#' regression is recorded as `alpha`, and the exponentiated intercept as `A`.
#'
#' Using `alpha`, `A`, `beta`, and `B`, estimates for `k21`, `kelim`, and `k12` can be made:
#'
#' - `k21 <- (A * beta + B*alpha) / (A + B)`
#' - `kelim <- (alpha * beta) / k21`
#' - `k12 <- alpha + beta - k21 - kel`
#'
#' Finally, the "absorption" phase data are predicted using the "late" linear
#' regression, and the residuals are calculated. Then, the "absorption"-phase
#' residuals are predicted using the "early"-residuals linear regression, and
#' the residuals of the residuals are calculated. A third linear regression is
#' performed on these "absorption" phase residuals-of-residuals. The negative
#' slope of this third regression provides an estimate for `kgutabs`.
#'
#' Using the estimates for `k21` and `kgutabs`, along with `A`, `alpha`, and
#' `beta`, an estimate can be made for `Fgutabs_V1`.
#'
#'  - `Fgutabs_V1 <- A * ( (kgutabs - alpha) * (beta - alpha))/( kgutabs * (k21 - alpha) )`
#'
#'  In case there is not enough data available in each phase to perform
#'  regression, the default or `httk`-derived starting values are used.
#'
#' #### Both PO and IV dosing data available
#'
#' The final starting values for `k21`, `kelim`, and `k12` are taken as the
#' average of the IV-based estimates and the PO-based estimates (with any NAs
#' removed).
#'
#' The estimates for `V1` (from the IV data) and `Fgutabs_V1` (from the PO data)
#' are multiplied to produce an estimate for `Fgutabs`.
#'
#'
#'
#' @param par_DF Optional: A data.frame as produced by [get_opt_params()], with
#'   a character variable `param_name` containing parameter names, and a logical
#'   variable `use_param` containing `TRUE` if parameter is to be fitted, and
#'   `FALSE` if parameter is not to be fitted. Optionally, may also contain
#'   numeric variables `lower_bound` and `upper_bound`, as produced by
#'   [get_lower_bounds()] and [get_upper_bounds()], respectively. Default NULL,
#'   in which case `par_DF` will be determined by calling [get_opt_params()]. If
#'   `par_DF` does not contain variables `lower_bound` and `upper_bound`, these
#'   will be added by calling [get_lower_bounds()] and [get_upper_bounds()].
#' @param model The name of the model whose parameters are to be estimated.
#' @param fitdata A `data.frame` or something that can be coerced to a
#'   `data.frame`: the concentration-time-dose data to be used for fitting. Must
#'   have harmonized variable names as produced by [rename_columns()].
#' @param starts_default A `data.frame` with variables `param_name` and
#'   `start_value`, defining default parameter starting values. For default
#'   value, see Details. Must include at least all variables relevant to the
#'   specified model with the dosing routes available in `fitdata` -- i.e.,
#'   those that would have `use_param == TRUE` in the `data.frame` returned by
#'   [get_opt_params()] with the specified `model` and `fitdata`.
#' @param start_from_httk Character vector: The names of one or more parameters
#'   whose starting values are to be determined using the function
#'   [httk::parameterize_1compartment]. Default value is `"all"`, which means
#'   that starting values for all parameters that *can* be estimated from
#'   [httk::parameterize_1compartment], *will* be estimated from
#'   [httk::parameterize_1compartment], so long as the `httk`-derived parameter
#'   values are within the lower and upper bounds defined by `par_DF$lower` and
#'   `par_DF$upper`. To suppress estimation of starting parameter values from
#'   `httk` data and just use the starting values in `starts_default`, set
#'   `start_from_httk` to NULL, NA, or a zero-length character vector.
#' @param start_from_data Character vector: The names of one or more parameters
#'   whose starting values are to be estimated from the data. Default value is
#'   `"all"`, which means that starting values for all parameters will be
#'   estimated from the data, so long as the `httk`-derived parameter
#'   values are within the lower and upper bounds defined by `par_DF$lower` and
#'   `par_DF$upper`. To suppress estimation of starting values from data, set
#'   `start_from_data` to NULL, NA, or a zero-length character vector. If both
#'   `start_from_httk` and `start_from_data` are `TRUE`, then data-derived
#'   starting values will override `httk` starting values.
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
                                         1), #sigma
                         start_value_msg = "Default"
                       ),
                       start_from_httk = "all",
                       start_from_data = "all",
                       suppress.messages = FALSE
                       ){

  fitdata <- as.data.frame(fitdata)

  if(all(is.null(start_from_data) |
     is.na(start_from_data) |
     length(start_from_data)==0 |
     !nzchar(start_from_data))){
    start_from_data <- NULL
  }

  if(all(tolower(start_from_data) %in% "all")){
    start_from_data <- par_DF$param_name
  }

  if(all(is.null(start_from_httk) |
         is.na(start_from_httk) |
         length(start_from_httk)==0 |
         !nzchar(start_from_httk))){
    start_from_httk <- NULL
  }

  if(all(tolower(start_from_httk) %in% "all")){
    start_from_httk <- par_DF$param_name
  }

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

    #the abovemerge won't get reference-specific sigmas: handle them
    for(this_sigma in grep(x = par_DF$param_name,
                           pattern= "sigma",
                           value = TRUE)){
      par_DF <- assign_start(param_name = this_sigma,
                             param_value = starts_default[
                               starts_default$param_name %in% "sigma",
                               "start_value"],
                             msg = "Default starting value",
                             par_DF = par_DF,
                             start_from = start_from_data)
    } #end for loop over reference-specific sigmas

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
        #set V1 to the httk-predicted Vdist if necessary
        #create a new "V1" parameter in httk_params for this purpose
        httk_params[["V1"]] <- httk_params[["Vdist"]]
      #replace the parameters in par_DF whose names match those in httk_params,
      #and are named in start_from_Httk
      #find matching names
      replace_names <- intersect(names(httk_params),
                                       intersect(par_DF$param_name,
                                                 start_from_httk))
      #update par_DF for any of these parameters
      for (this_param in replace_names){
        if(this_param %in% par_DF$param_name){
          par_DF <- assign_start(param_name = this_param,
                                 param_value = httk_params[[this_param]],
                                 msg = paste0("httk::parameterize_1comp(",
                                              "DTXSID = ", this.dtxsid, ", ",
                                              "default.to.human = TRUE, ",
                                              "suppress.messages = TRUE",
                                              "species = ", httk_species,
                                              ")"),
                                 start_from = start_from_httk,
                                 par_DF <- par_DF)
        }
      }


    } #end if(this.dtxsid %in% httk::get_cheminfo(info="dtxsid"))

    #Try roughly estimating parameters from data

    #Drop control points & non-detects
    tmpdat <- subset(fitdata,
                     Dose > 0 &
                       is.finite(Value))

    #normalize concentration by dose
    tmpdat$ValueDose <- tmpdat$Value/tmpdat$Dose
    #log transform Value/Dose
    tmpdat$logValueDose <- log(tmpdat$ValueDose)

    #if model is flat, take A to be the mean log concentration/dose, on natural scale
    if(model %in% "flat"){
      A <- exp(mean(tmpdat$logValueDose, na.rm = TRUE))
      par_DF <- assign_start(param_name = "A",
                             param_value = A,
                             msg = "Median concentration/dose",
                             start_from = start_from_data,
                             par_DF = par_DF)
    }else{


    #Split into IV and oral datasets
    iv_data <- tmpdat[tmpdat$Route %in% "iv", ] #will be empty if no IV data
    po_data <- tmpdat[tmpdat$Route %in% "po", ] #will be empty if no PO data

    has_iv <- any(tmpdat$Route %in% "iv")
    has_po <- any(tmpdat$Route %in% "po")

    #####################
    # 1-compartment model
    #####################
    if(model %in% "1compartment"){
      if(has_iv %in% TRUE){
        #--------------------------------------------------
        #get kelim and Vdist from linear regression
        #see https://www.boomer.org/c/p4/c04/c0406.php
        #--------------------------------------------------
        if(length(unique(iv_data$Time))>=2){
          if(nrow(iv_data)>2){
            #if there is enough data to do linear regression, get Vdist and kelim from a
            #linear regression of log (Value/Dose) vs. Time
            lm_iv <- lm(logValueDose ~ Time,
                        data = iv_data)
            Vdist_iv <- 1/exp(coef(lm_iv)[1])
            kelim_iv <- -coef(lm_iv)[2]
          }else if(nrow(iv_data)==2){
            #if there are exactly 2 data points,
            iv_data <- iv_data[order(iv_data$Time), ]
            #draw a line connecting them
            kelim_iv <- - diff(iv_data$logValueDose)/diff(iv_data$Time)
            #extrapolate back to time 0 to get Vdist
            log_A_Dose <- iv_data[1, "logValueDose"] + kelim_iv * iv_data[1, "Time"]
            Vdist_iv <- 1/exp(log_A_Dose)
          }
        }else{
          #if there is only 1 IV data point, or none
          Vdist_iv <- NA_real_
          kelim_iv <- NA_real_
        }

        par_DF <- assign_start(param_name = "Vdist",
                               param_value = Vdist_iv,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Linear regression of IV data")
        if(!is.finite(Vdist_iv) |
           Vdist_iv < par_DF["Vdist", "lower_bound"]){
          #if trying to fit Vdist failed,
          #then try to get Vdist as 1/median Value/Dose at min time
          Vdist_iv <- 1/median(iv_data[iv_data$Time == min(iv_data$Time),
                                       "ValueDose"])
          par_DF <- assign_start(param_name = "Vdist",
                                 param_value = Vdist_iv,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "1/(Median Value/Dose) at min time")
        }

        #if trying to fit kelim failed,
        #then draw a line between the last IV time point and log(1/Vdist)
        if(!is.finite(kelim_iv) | kelim_iv < par_DF["kelim", "lower_bound"]){
          iv_data <- iv_data[order(iv_data$Time, decreasing = TRUE), ]
          kelim_iv <- -(iv_data[1, "logValueDose"] - log(1/Vdist_iv))/
            (iv_data[1, "Time"] - min(iv_data$Time))

          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kelim_iv,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Slope of line between log Value/Dose for min and max time points in IV data")
        }

      }else{
        #if no iv data, just keep NA as estimates
        kelim_iv <- NA_real_
        Vdist_iv <- NA_real_
      }

      if(has_po %in% TRUE){
        #--------------------------------------------------
        #get kelim and Fgutabs/Vdist from linear regression
        #see https://www.boomer.org/c/p4/c09/c0901.php
        #--------------------------------------------------
        #split data into early and late parts,
        #using the time of peak concentration as the dividing line

        #to get time of peak concentration:
        #get average log conc/dose by time point
        Cavg <- tapply(po_data$logValueDose,
                       po_data$Time,
                       FUN = function(x) mean(x, na.rm = TRUE))
        #Cavg will be sorted by ascending time
        #get corresponding time points for the averages
        Cavg_times <- as.numeric(names(Cavg))

        Cpeak <- max(Cavg, na.rm = TRUE) #peak average log conc/dose
        tpeak <- Cavg_times[Cavg == Cpeak] #time of peak log conc/dose

        #if tpeak is the first or last time, bump it to the second or second to last
        if(all(po_data$Time >= tpeak)){ #if it is the first time
          tpeak <- Cavg_times[2]
          Cpeak <- Cavg[2]
        }else if(all(po_data$Time <= tpeak)){ #if it is the last time
          tpeak <- Cavg_times[length(Cavg_times)-1]
          Cpeak <- Cavg[length(Cavg)-1]
        }

        po_early <- subset(po_data,
                           Time <= tpeak)
        po_late <- subset(po_data,
                          Time >= tpeak) #yes, include tpeak in both

        if(length(unique(po_late$Time)) >= 2){
          if(nrow(po_late) > 2){ #if enough late data to do regression
            #linear regression of late data
            lm_late <- lm(logValueDose ~ Time, data = po_late)
            kelim_po <- -coef(lm_late)[2]
            A_Dose_po <- exp(coef(lm_late)[1]) #intercept -- use later to get Fgutabs/Vdist

          }else{ #if only 2 points in late data
            po_late <- po_late[order(po_late$Time), ]
            #draw a line between points
            kelim_po <- -diff(po_late$logValueDose)/diff(po_late$Time)
            #extrapolate back to time 0
            A_Dose_po <- exp(po_late[1, "logValueDose"] + kelim_po * po_late[1, "Time"])
          }
        }else{ #if only one point in late data, or none
          A_Dose_po <- NA_real_
          kelim_po <- NA_real_
        }

        #if regression failed for kelim, then assume max time = 5 * elimination half-life
        if(!is.finite(kelim_po) | kelim_po <= par_DF["kelim", "lower_bound"]){
          thalf_elim <- max(po_data$Time)/5
          kelim_po <- log(2)/thalf_elim
          #extrapolate back to time 0
          A_Dose_po <- exp(po_late[1, "logValueDose"] + kelim_po * po_late[1, "Time"])
        }

        if(has_iv %in% TRUE){  #if we also have iv data, take average kelim value
          par_DF <- assign_start(param_name = "kelim",
                                 param_value =  mean(c(kelim_iv, kelim_po), na.rm = TRUE),
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Avg of lin reg on IV data and method of residuals on PO data")
        }else{
          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kelim_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data")
        }

        #get residuals
        #predicted values:
        pred_late <- exp(log(A_Dose_po) - kelim_po * po_early$Time)
        #residuals: late predictions - early observations
        resid_early <- pred_late - po_early$ValueDose
        #log-transform residuals
        suppressWarnings(po_early$logresid <- log(resid_early))
        #drop any NA log residuals
        po_early <- subset(po_early, is.finite(logresid))

        if(length(unique(po_early$Time)) >= 2){
          if(nrow(po_early) > 2){
            #regress log residuals on Time
            lm_resid <- lm(logresid ~ Time,
                           data = po_early)
            kgutabs_po <- -coef(lm_resid)[2]
          }else if(nrow(po_early)==2){ #if only 2 points, draw a line through them
            po_early <- po_early[order(po_early$Time), ]
            kgutabs_po <- -(diff(po_early$logresid)/diff(po_early$Time))
          }
        }else{ #if 1 or 0 non-NA residuals
          kgutabs_po <- NA_real_
        }

        #update par_DF
        #kgutabs
        par_DF <- assign_start(param_name = "kgutabs",
                               param_value = kgutabs_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on PO data")

        #if method of residuals failed, then try just assuming that tpeak = 5 *
        #absorption half-life -- see https://www.boomer.org/c/p4/c09/c0904.php
        if(!is.finite(kgutabs_po) | kgutabs_po < par_DF["kgutabs", "lower_bound"]){
          thalf_abs <- tpeak/5
          kgutabs_po <- log(2)/thalf_abs
          #update par_DF
          #kgutabs
          par_DF <- assign_start(param_name = "kgutabs",
                                 param_value = kgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of inspection on PO data (assume tpeak = 5 * absorp half-life)")

        }



        #use intercept to get Fgutabs_Vdist
        Fgutabs_Vdist_po <- A_Dose_po * (kgutabs_po - kelim_po) / kgutabs_po

        #Fgutabs_Vdist
        par_DF <- assign_start(param_name = "Fgutabs_Vdist",
                               param_value = Fgutabs_Vdist_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on PO data")

        if(has_iv %in% TRUE){  #if we also have iv data
          #then we can estimate Fgutabs separately
          Vdist <- par_DF["Vdist", "start_value"]
          Fgutabs_po <- Fgutabs_Vdist_po * Vdist

          #update par_DF
          #Fgutabs
          par_DF <- assign_start(param_name = "Fgutabs",
                                 param_value = Fgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data to get Fgutabs/Vdist, combined with estimate of Vdist")


        }
      } #end if(has_po %in% TRUE)
    } #end if(model %in% "1compartment")

    #####################
    # 2-compartment model
    #####################
    if(model %in% "2compartment"){
      if(has_iv %in% TRUE & length(unique(iv_data$Time)) > 2){
        #--------------------------------------------------
        #get alpha, beta, A, and B from method of residuals
        #see https://www.boomer.org/c/p4/c19/c1903.php
        #--------------------------------------------------

        #split IV data into early and late parts
        #find the dividing line between early and late as the "elbow" point
        elbow <- tryCatch(akmedoids::elbow_point(iv_data$Time,
                                             iv_data$logValueDose)[c("x", "y")],
                               error = function(err){
                                 #if akmedoids::elbow_point() does not work,
                                 #fallback to middle time (naive)
                                 elbow_x <- (max(iv_data$Time) - min(iv_data$Time))/2
                                 elbow_y <- approx(x = iv_data$Time,
                                                   y = iv_data$logValueDose,
                                                   xout = elbow_x)
                                 return(list(x = elbow_x,
                                               y = elbow_y))
                               })

        if(!is.finite(elbow$x) |
           elbow$x < 0){
          elbow_x <- (max(iv_data$Time) - min(iv_data$Time))/2
          elbow_y <- approx(x = iv_data$Time,
                            y = iv_data$logValueDose,
                            xout = elbow_x)
          elbow <- list(x = elbow_x,
                      y = elbow_y)
        }

        elbow_time <- elbow$x

          #if elbow time does not allow 2 time points in both early and late phases,
          #then move it so that it does
          if(sum(iv_data$Time >= elbow_time)<2){
            elbow_x <- sort(unique(iv_data$Time))[2]
            elbow_y <- approx(x = iv_data$Time,
                              y = iv_data$logValueDose,
                              xout = elbow_x)
            elbow <- list(x = elbow_x,
                          y = elbow_y)
            elbow_time <- elbow$x
          }else if(sum(iv_data$Time <= elbow_time)<2){
            elbow_x <- sort(unique(iv_data$Time),
                            decreasing = TRUE)[2]
            elbow_y <- approx(x = iv_data$Time,
                              y = iv_data$logValueDose,
                              xout = elbow_x)
            elbow <- list(x = elbow_x,
                          y = elbow_y)
            elbow_time <- elbow$x
          }

        #split IV data into early and late phases at elbow_time
        iv_early <- subset(iv_data,
                           Time <= elbow_time)
        iv_late <- subset(iv_data,
                          Time >= elbow_time)

        #--------------------------------------------------
        # Regression on late-phase data
        #-------------------------------------------------

        if(length(unique(iv_late$Time))>=2){
          if(nrow(iv_late) > 2){
            #if enough late data to regress
            #linear regression of late data gives beta, B/Dose
            lm_late <- lm(logValueDose ~ Time, data = iv_late)
            beta_iv <- -coef(lm_late)[2]
            B_Dose_iv <- exp(coef(lm_late)[1])
          }else{
            #if we have exactly 2 time points, draw a line through them
            iv_late <- iv_late[order(iv_late$Time), ]
            beta_iv <- -diff(iv_late$logValueDose)/diff(iv_late$Time)
            #intercept: extrapolate back to time zero
            B_Dose_iv <- exp(iv_late[1, "logValueDose"] + beta_iv * iv_late[1, "Time"])
          }
        }else{ #if we don't have at least two late time points, can't do anything
          B_Dose_iv <- NA_real_
          beta_iv <- NA_real_
        }

        #if regression failed for beta, assume max time is 5 * terminal half-life
        if(!is.finite(beta_iv) | beta_iv <= 0){
          #use max time for all IV data, whether detect or no
          iv_all <- subset(fitdata, Route %in% "iv")
          thalf_beta <- max(iv_all$Time)/5
          beta_iv <- log(2)/thalf_beta
          #extrapolate back from elbow point to zero
          B_Dose_iv <- exp(elbow$y + beta_iv * elbow$x)
        }

        #---------------------------------------------------
        # Regression on early-phase RESIDUALS
        #---------------------------------------------------

        #calc residuals during early phase observed - (predicted from late-phase
        #regression) (because late-phase predictions should be less than
        #early-phase observations)
        pred_late <- exp(log(B_Dose_iv) + -beta_iv * iv_early$Time)
        resid_early <- iv_early$ValueDose - pred_late
        suppressWarnings(iv_early$logresid <-  log(resid_early))
        #drop any NA residuals
        iv_early <- subset(iv_early, is.finite(logresid))

        #regress log residuals on Time
        if(length(unique(iv_early$Time))>=2){
        if(nrow(iv_early) > 2){ #if enough early-phase data for regression
          lm_resid <- lm(logresid ~ Time,
                         data = iv_early)
        alpha_iv <- -coef(lm_resid)[2]
        A_Dose_iv <- exp(coef(lm_resid)[1])
        }else if(nrow(iv_early) == 2){ #if only 2 early-phase points
            #get slope & intercept by drawing line through 2 points
            iv_early <- iv_early[order(iv_early$Time), ]
            alpha_iv <- -(diff(iv_early$logresid)/diff(iv_early$Time))
            #extrapolate back from beginning of early-phase to time = zero
            A_Dose_iv <- exp(iv_early[1, "logValueDose"] + alpha_iv * iv_early[1, "Time"])
        }
      }else{ #if one or zero early-phase time points with finite log residuals
          A_Dose_iv <- NA_real_
          alpha_iv <- NA_real_
        }

        #if method of residuals failed for alpha,
        #then try assuming that elbow point = 5 * early-phase half life
        if(!is.finite(alpha_iv) | alpha_iv <= 0){
          thalf_alpha_iv <- elbow_time/5
          alpha_iv <- log(2)/thalf_alpha_iv
          #use iv_early data before NA residuals were dropped
          #extrapolate back from elbow point to time = zero
          A_Dose_iv <- exp(elbow$y + alpha_iv * elbow$x)
        }

        k21_iv <- (A_Dose_iv * beta_iv + B_Dose_iv*alpha_iv)/(A_Dose_iv + B_Dose_iv)
        kel_iv <- (alpha_iv * beta_iv)/k21_iv
        k12_iv <- alpha_iv + beta_iv - k21_iv - kel_iv
        V1_iv <- 1/(A_Dose_iv + B_Dose_iv)
        #to show V1 relation:
        #see https://www.boomer.org/c/p4/c19/c1902.php for A. B equations;
        #write A+B and solve for V1

        #update par_DF
        par_DF <- assign_start(param_name = "V1",
                               param_value = V1_iv,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on IV data")

        #k21
        par_DF <- assign_start(param_name = "k21",
                               param_value = k21_iv,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on IV data")
        #kelim
        par_DF <- assign_start(param_name = "kelim",
                               param_value = kel_iv,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on IV data")
        #k12
        par_DF <- assign_start(param_name = "k12",
                               param_value = k12_iv,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on IV data")
      } #end  if(has_iv %in% TRUE)

      #----------------------------------------------------------------
      # Oral data -- 2 compartment model
      # See https://www.boomer.org/c/p4/c19/c1907.php
      # Use method of residuals *twice*
      #----------------------------------------------------------------

      if(has_po %in% TRUE){
        #Split the data into three parts:
        #Absorption phase, early phase, and late phase.

        #Absorption phase ends when peak concentration is achieved:
        #To get time of peak concentration:
        #get average log conc/dose by time point
        Cavg <- tapply(po_data$logValueDose,
                       po_data$Time,
                       FUN = function(x) mean(x, na.rm = TRUE))
        #get corresponding time points for the averages
        Cavg_times <- as.numeric(names(Cavg))

        Cpeak <- max(Cavg, na.rm = TRUE) #peak average conc (in units of log Conc/Dose)
        tpeak <- Cavg_times[Cavg == Cpeak] #time of peak conc

        #split data into absorption and non-absorption phases
        po_abs <- subset(po_data, Time <= tpeak)
        po_nonabs <- subset(po_data, Time > tpeak)

        #further split the non-absorption phase into early and late parts
        #find the dividing line between early and late as the "elbow" point
        elbow <- tryCatch(akmedoids::elbow_point(po_nonabs$Time,
                                           po_nonabs$logValueDose)[c("x", "y")],
                               error = function(err){
                                 #if akmedoids::elbow_point() does not work
                                 #fall back to middle time (naive)
                                 elbow_x <- (max(po_nonabs$Time) - min(po_nonabs$Time))/2
                                 elbow_y <- approx(x = po_nonabs$Time,
                                                   y = po_nonabs$logValueDose,
                                                   xout = elbow_x)
                                return(list(x = elbow_x,
                                      y = elbow_y))
                               })



        #in case we get NA or infinite or negative elbow time,
        #fallback to middle time
        if(!is.finite(elbow$x) |
           elbow$x < 0){
          elbow_x <- (max(po_nonabs$Time) - min(po_nonabs$Time))/2
          elbow_y <- approx(x = po_nonabs$Time,
                            y = po_nonabs$logValueDose,
                            xout = elbow_x)
          elbow <- list(x = elbow_x,
                        y = elbow_y)
        }

        elbow_time <- elbow$x

        if(length(unique(po_nonabs$Time))>1){
        #if elbow time does not allow 2 points in both early and late phases,
          #then move it so that it does
        if(sum(po_nonabs$Time >= elbow_time)<2){
          elbow_x <- sort(unique(po_nonabs$Time))[2]
          elbow_y <- approx(x = po_nonabs$Time,
                            y = po_nonabs$logValueDose,
                            xout = elbow_x)
          elbow <- list(x = elbow_x,
                        y = elbow_y)
          elbow_time <- elbow$x
        }else if(sum(po_nonabs$Time <= elbow_time)<2){
          elbow_x <- sort(unique(po_nonabs$Time),
                             decreasing = TRUE)[2]
          elbow_y <- approx(x = po_nonabs$Time,
                            y = po_nonabs$logValueDose,
                            xout = elbow_x)
          elbow <- list(x = elbow_x,
                        y = elbow_y)
          elbow_time <- elbow$x
        }
        }

        po_early <- subset(po_nonabs, Time <= elbow_time)
        po_late <- subset(po_nonabs, Time >= elbow_time)

        if(length(unique(po_late$Time))>=2){
          if(nrow(po_late) > 2){
            #linear regression of late data gives beta, B/Dose
            lm_late <- lm(logValueDose ~ Time, data = po_late)
            beta_po <- -coef(lm_late)[2]
            B_Dose_po <- exp(coef(lm_late)[1])
          }else if(nrow(po_late) == 2){
            #if only 2 timepoints of late data, draw a line through them
            po_late <- po_late[order(po_late$Time), ]
            beta_po <- -(diff(po_late$logValueDose)/diff(po_late$Time))
            #intercept -- extrapolate back to time 0
            B_Dose_po <- exp(po_late[1, "logValueDose"] + beta_po * po_late[1, "Time"])
          }
        }else{ #if one or zero points
          B_Dose_po <- NA_real_
          beta_po <- NA_real_
        }

        #if regression fails for beta_po,
        #assume last time point is 5 times terminal half-life
        if(!is.finite(beta_po) |
           beta_po < 0){
          #use max time for all po data, whether detect or nondetect
          po_all <- subset(fitdata, Route == "po")
          thalf_beta <- max(po_all$Time)/5
          beta_po <- log(2)/thalf_beta
          #intercept: extrapolate back to time = 0 from elbow point
          B_Dose_po <- exp(elbow$y + beta_po * elbow_time)
        }

        #residuals for early data
        #obs should be greater than predicted
        #calc residuals during early phase observed - (predicted from late-phase
        #regression) (because late-phase predictions should be less than
        #early-phase observations)
        pred_late_po <- exp(log(B_Dose_po) + -beta_po * po_early$Time)
        resid_early_po <- po_early$ValueDose - pred_late_po
        suppressWarnings(po_early$logresid <-  log(resid_early_po))
      #drop any NA log residuals
        po_early <- subset(po_early,
                           is.finite(logresid))

        if(length(unique(po_early$Time))>=2){
          if(nrow(po_early) > 2){
            #linear regression of residuals gives alpha, A/Dose
            lm_resid_early <- lm(logresid ~ Time,
                                 data = po_early)
            alpha_po <- -coef(lm_resid_early)[2]
            A_Dose_po <- exp(coef(lm_resid_early)[1])
          }else if(nrow(po_early) == 2){ #if 2 points

            po_early <- po_early[order(po_early$Time), ]
            alpha_po <- -(diff(po_early$logresid)/diff(po_early$Time))
            A_Dose_po <- exp(po_early[1, "logresid"])
          }
        }else{ #if one or zero time points
          alpha_po <- NA_real_
          A_Dose_po <- NA_real_
        }

        #if method of residuals failed for alpha,
        #then try assuming that elbow point = 5 * early-phase half life
        if(!is.finite(alpha_po) | alpha_po <= 0){
          thalf_alpha <- elbow_time/5
          alpha_po <- log(2)/thalf_alpha
          #extrapolate from elbow point back to time = 0
            A_Dose_po <- exp(elbow$y + alpha_po * elbow_time)
        }

        #residuals for absorption phase
        #predictions from late-phase regression:
        late_pred <- exp(log(B_Dose_po) + -beta_po * po_abs$Time)
        #predictions from early-phase residuals regression
        early_pred <- exp(log(A_Dose_po) + -alpha_po * po_abs$Time)

        #predictions should be > observations now, so reverse sign to calc resids
        po_abs$resid <- (late_pred + early_pred) - po_abs$ValueDose
        suppressWarnings(po_abs$logresid <- log(po_abs$resid))

        #drop NA residuals
        po_abs <- subset(po_abs,
                         is.finite(logresid))

        if(length(unique(po_abs$Time))>=2){
          if(nrow(po_abs) > 2){
            #linear regression of residuals gives kgutabs
            lm_resid_abs <- lm(logresid ~ Time,
                               data = po_abs)

            kgutabs_po <- -coef(lm_resid_abs)[2]
          }else{
            po_abs <- po_abs[order(po_abs$Time), ]
            kgutabs_po <- -(diff(po_abs$logresid)/diff(po_abs$Time))
          }
        }else{ #if zero or 1 time points
          kgutabs_po <- NA_real_
        }

        par_DF <- assign_start(param_name = "kgutabs",
                               param_value = kgutabs_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on PO data")

        #if method of residuals failed, then try just assuming that tpeak = 5 *
        #absorption half-life -- see https://www.boomer.org/c/p4/c09/c0904.php
        if(!is.finite(kgutabs_po) | kgutabs_po <= 0){
          thalf_abs <- tpeak/5
          kgutabs_po <- log(2)/thalf_abs
          par_DF <- assign_start(param_name = "kgutabs",
                                 param_value = kgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Assuming tpeak = 5 * absorption half-life")
        }



        #get k21, k12, kel from A, B, alpha, beta
        #source: https://www.boomer.org/c/p4/c19/c1903.php
        k21_po <- (A_Dose_po * beta_po + B_Dose_po*alpha_po) / (A_Dose_po + B_Dose_po)
        kel_po <- (alpha_po * beta_po) / k21_po
        k12_po <- alpha_po + beta_po - k21_po - kel_po

        #Solve A_Dose for Fgutabs/V1
        #See cp_2comp()
        Fgutabs_V1_po <- A_Dose_po *
          ( (kgutabs_po - alpha_po) * (beta_po - alpha_po))/
          ( kgutabs_po * (k21_po - alpha_po) )
        #update par_DF with Fgutabs_V1
        par_DF <- assign_start(param_name = "Fgutabs_V1",
                               param_value = Fgutabs_V1_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on PO data")

        if(has_iv %in% TRUE){ #if both po and iv data
          #then we already have an estimate of V1, and can estimate Fgutabs by itself
          Fgutabs_po <- Fgutabs_V1_po * V1_iv
          #update par_DF with Fgutabs
          par_DF <- assign_start(param_name = "Fgutabs",
                                 param_value = Fgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data with V1 from method of residuals on IV data")
          #k21, kelim, k12:
          #update par_DF using average of PO and IV estimates
          k21_avg <- mean(c(k21_iv,
                            k21_po),
                          na.rm = TRUE)
          par_DF <- assign_start(param_name = "k21",
                                 param_value = k21_avg,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Avg of method of residuals on PO and IV data")
          kel_avg <- mean(c(kel_iv,
                            kel_po),
                          na.rm = TRUE)
          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kel_avg,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Avg of method of residuals on PO and IV data")
          k12_avg <- mean(c(k12_iv,
                            k12_po),
                          na.rm = TRUE)
          par_DF <- assign_start(param_name = "k12",
                                 param_value = k12_avg,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Avg of method of residuals on PO and IV data")
        }else{ #if no IV data, only PO data
          #k21, kelim, k12:
          #update par_DF using PO estimates
          par_DF <- assign_start(param_name = "k21",
                                 param_value = k21_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data")

          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kel_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data")

          par_DF <- assign_start(param_name = "k12",
                                 param_value = k12_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data")
        } #end if(has_iv %in% TRUE)/else
      } #end if(has_po %in% TRUE)
    } #end if(model %in% "2compartment")
    }#end if model %in% "flat"/else

    #sigma:
  #Try evaluating model with the starting values for parameters
  #Then get SD of residuals
  params <- as.list(par_DF[par_DF$use_param %in% TRUE, "start_value"])
  names(params) <- par_DF[par_DF$use_param %in% TRUE, "param_name"]
  if(model %in% "1compartment"){
    modelf <- "cp_1comp"
  }else if(model %in% "2compartment"){
    modelf <- "cp_2comp"
  }else if(model %in% "flat"){
    modelf <- "cp_flat"
  }

  pred <- do.call(modelf,
                  list(params = params,
                       time = fitdata$Time,
                       dose = fitdata$Dose,
                       iv.dose = fitdata$Route %in% "iv"))

  resid <- pred - fitdata$Value
  resid[!is.finite(resid)] <- NA_real_

  for(this_sigma in grep(x = par_DF$param_name,
                         pattern= "sigma",
                         value = TRUE)){
    if(this_sigma == "sigma"){ #if only one sigma (one reference, or pooled)
      tmp_sigma <- sd(resid, na.rm = TRUE)
    }else{
      #if multiple sigmas for multiple references
      #get the reference
      refid <- gsub(x = this_sigma,
                    pattern = "sigma_ref_",
                    replacement = "")
      resid_ref <- resid[fitdata$Reference %in% refid]
      tmp_sigma <- sd(resid_ref, na.rm = TRUE)
    }

    par_DF <- assign_start(param_name = this_sigma,
                           param_value = tmp_sigma,
                           msg = "SD of resid for model with starting values",
                           par_DF = par_DF,
                           start_from = start_from_data)
  } #end for loop over sigmas

  #For anything with use_param == FALSE, set the starting value to NA
    #because no starting value will be used
  par_DF[!(par_DF$use_param %in% TRUE),
         c("start_value",
           "start_value_msg")] <- list(NA_real_,
                                    "use_param is not TRUE")

  return(par_DF)
}

#' Assign starting value
#'
#' Utility function to assign starting value for a parameter, if within
#' specified bounds
#'
#' @param param_name Character: The name of the parameter whose starting value
#'   should be checked and updated
#' @param param_value Numeric: The proposed starting value for the parameter
#' @param msg Character: The proposed message to use if the proposed starting
#'   value is acepted
#' @param par_DF Data frame: A `data.frame` of parameters and bounds, as
#'   produced by [get_opt_params], [get_lower_bounds], and [get_upper_bounds].
#'   Must have row names corresponding to parameter names, and must have
#'   variables `lower_bound` and `upper_bound`.
#' @param start_from Character: A vector of parameter names for which to
#'   assign starting values from data. May be "all" to assign all from data, or
#'   NULL, NA, character(0), or `''` to assign none from data.
#'
#' @return If the proposed starting value is within the bounds for the specified
#'   parameter, this function returns `par_DF` updated with variables
#'   `start_value` and `start_value_msg` in the row for the specified parameter,
#'   now containing `param_value` and `msg`, respectively. If those variables did
#'   not already exist in `par_DF`, they will be added. If the proposed starting
#'   value is *not* within bounds, this function returns `par_DF` unchanged.
#' @author Caroline Ring

assign_start <- function(param_name,
                         param_value,
                         msg,
                         par_DF,
                         start_from){
  #get lower bounds
  lower <- par_DF[par_DF$param_name %in% param_name, "lower_bound"]
  upper <- par_DF[par_DF$param_name %in% param_name, "upper_bound"]

  if(all(is.null(start_from) |
         is.na(start_from) |
         length(start_from)==0 |
         !nzchar(start_from))){
    start_from <- NULL
  }

  if(all(tolower(start_from) %in% "all")){
    start_from <- par_DF$param_name
  }

  if(!is.null(param_value)){
  if(is.finite(param_value) &
    (param_value >= lower) %in% TRUE &
     (param_value <= upper) %in% TRUE &
     (param_name %in% start_from) %in% TRUE){     #if within bounds and finite, update par_DF
    par_DF[par_DF$param_name %in% param_name, c("start_value",
                         "start_value_msg")] <- list(param_value,
                                                     msg)
  } #if outside bounds or not finite, return par_DF unchanged
  }

  return(par_DF)
}
