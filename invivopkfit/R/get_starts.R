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
#' | sigma          | 100         | Default         |
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
#' phases, naively using the median value of Time as the dividing line.
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
#' divided into "early" and "late" parts at the median value of Time.
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
#' ### PO dosing data available
#'
#' `kelim`, `k12`, `k21` `kgutabs`, and `Fgutabs_V1` (the ratio of `Fgutabs` to
#' `V1`) are estimated using the method of residuals on `log(Value/Dose)` vs.
#' `Time` from the PO data. Specifically, the data are divided into three parts:
#' "absorption" phase, "early" phase, and "late" phase, naively using the 25th
#' percentile of Time as the boundary between "absorption" and "early" phases,
#' and the 75th percentile of Time as the boundary between "early" and "late"
#' phases.
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
#' @param par_DF Optional: A data.frame as produced by [get_opt_params], with a
#'   character variable [param_name] containing parameter names, and a logical
#'   variable [use_param] containing TRUE if parameter is to be fitted, and
#'   FALSE if parameter is not to be fitted. Optionally, may also contain
#'   numeric variables `lower_bound` and `upper_bound`, as produced by
#'   [get_lower_bounds] and [get_upper_bounds], respectively. Default NULL, in
#'   which case `par_DF` will be determined by calling [get_opt_params]. If
#'   `par_DF` does not contain `lower_bound` and `upper_bound`, these will be
#'   added by calling [get_lower_bounds] and [get_upper_bounds].
#' @param model The name of the model whose parameters are to be estimated.
#' @param fitdata A `data.frame` or something that can be coerced to a
#'   `data.frame`: the concentration-time-dose data to be used for fitting. Must
#'   have harmonized variable names as produced by [rename_columns].
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
                                         100), #sigma
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





    #if model is flat, take A to be the average concentration
    if(model %in% "flat"){
      A <- mean(fitdata$Value, na.rm = TRUE)
      par_DF <- assign_start(param_name = "A",
                             param_value = A,
                             msg = "Mean concentration",
                             start_from = start_from_data,
                             par_DF = par_DF)
    }

    #Split into IV and oral data
    tmpdat <- subset(fitdata, !is.na(Value) &
                       Dose > 0)

    iv_data <- tmpdat[tmpdat$Route %in% "iv", ] #will be empty if no IV data
    po_data <- tmpdat[tmpdat$Route %in% "po", ] #will be empty if no PO data

    has_iv <- any(tmpdat$Route %in% "iv")
    has_po <- any(tmpdat$Route %in% "po")

    #normalize concentration by dose
    iv_data$ValueDose <- iv_data$Value / iv_data$Dose
    po_data$ValueDose <- po_data$Value / po_data$Dose

    if(model %in% "1compartment"){
      if(has_iv %in% TRUE){
        #get Vdist, kelim from a linear regression of log (Value/Dose) vs. Time
        Vdist_iv <- exp(coef(lm(log(ValueDose) ~ Time,
                                data = iv_data))[1])
        kelim_iv <- -coef(lm(log(ValueDose) ~ Time,
                             data = iv_data))[2]

        #update par_DF, checking if estimates are within bounds
        par_DF <- assign_start(param_name = "Vdist",
                               param_value = Vdist_iv,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Linear regression of IV data")

        par_DF <- assign_start(param_name = "kelim",
                               param_value = kelim_iv,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Linear regression of IV data")
      }

      if(has_po %in% TRUE){
        #method of residuals:
        #split data into early and late parts
        po_early <- subset(po_data,
                           Time <= median(Time))
        po_late <- subset(po_data,
                          Time > median(Time))
        #linear regression of late data
        lm_late <- lm(log(ValueDose) ~ Time, data = po_late)
        #kelim is negative slope of late data
        kelim_po <--coef(lm_late)[2]
        A_Dose_po <- exp(coef(lm_late)[1]) #intercept -- use later to get Fgutabs/Vdist
        #log-scale residuals
        po_early$logresid <- log(po_early$ValueDose) -
          predict(lm_late, newdata = po_early)
        #regress log residuals on Time again
        lm_resid <- lm(logresid ~ Time,
                       data = po_early)
        kgutabs_po <- -coef(lm_resid)[2]
        #use intercept to get Fgutabs_Vdist
        Fgutabs_Vdist_po <- A_Dose_po * (kgutabs_po - kelim_po) / kgutabs_po

        #update par_DF
        #kgutabs
        par_DF <- assign_start(param_name = "kgutabs",
                               param_value = kgutabs_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on PO data")


        #Fgutabs_Vdist
        par_DF <- assign_start(param_name = "Fgutabs_Vdist",
                               param_value = Fgutabs_Vdist_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on PO data")


        if(has_iv %in% TRUE){  #if we also have iv data
          #then we can estimate Fgutabs separately
          Fgutabs_po <- Fgutabs_Vdist_po * Vdist_iv

          #update par_DF
          #Fgutabs
          par_DF <- assign_start(param_name = "Fgutabs",
                                 param_value = Fgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data")

          #kelim
          #average IV and PO estimates
          kelim_avg <- mean(kelim_iv,
                            kelim_po,
                            na.rm = TRUE)
          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kelim_avg,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Avg of lin reg on IV data and method of residuals on PO data")

        }else{ #if no IV data
          #update par_DF with kelim from PO data
          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kelim_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data")

        } #end if(has_iv %in% TRUE)/else
      } #end if(has_po %in% TRUE)
    } #end if(model %in% "1compartment")

    #####################
    # 2-compartment model
    #####################
    if(model %in% "2compartment"){
      if(has_iv %in% TRUE){
        #get alpha, beta, A, and B from method of residuals
        #see https://www.boomer.org/c/p4/c19/c1903.php
        #split IV data into early and late parts
        iv_early <- subset(iv_data,
                           Time <= median(Time))
        iv_late <- subset(iv_data,
                          Time > median(Time))
        #linear regression of late data gives beta, B/Dose
        lm_late <- lm(log(ValueDose) ~ Time, data = iv_late)
        beta_iv <- -coef(lm_late)[2]
        B_Dose_iv <- exp(coef(lm_late)[1])

        #log-scale residuals
        iv_early$logresid <- log(iv_early$ValueDose) -
          predict(lm_late, newdata = iv_early)
        #regress log residuals on Time again
        lm_resid <- lm(logresid ~ Time,
                       data = iv_early)
        alpha_iv <- -coef(lm_resid)[2]
        A_Dose_iv <- exp(coef(lm_resid)[1])

        k21_iv <- (A_Dose_iv * beta_iv + B_Dose_iv*alpha_iv)/(A_Dose_iv + B_Dose_iv)
        kel_iv <- (alpha_iv * beta_iv)/k21_iv
        k12_iv <- alpha_iv + beta_iv - k21_iv - kel_iv
        #solve A for V1, see https://www.boomer.org/c/p4/c19/c1902.php
        #(and remember we already normalized by Dose)
        V1_iv <- (alpha_iv - k21_iv)/(A_Dose_iv * (alpha_iv - beta_iv))

        #update par_DF
        #V1
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

      if(has_po %in% TRUE){
        #run method of residuals twice to get kgutabs, alpha, beta, A, and B
        #this time, split the data into three parts:
        #absorption phase, early phase, and late phase.
        #Naively, just split time into tertiles.
        po_abs <- subset(po_data, Time < quantile(Time, 0.25))
        po_early <- subset(po_data, Time >= quantile(Time, 0.25) &
                             Time <quantile(Time, 0.75))
        po_late <- subset(po_data, Time >= quantile(Time, 0.75))

        #linear regression of late data gives beta, B/Dose
        lm_late <- lm(log(ValueDose) ~ Time, data = po_late)
        beta_po <- -coef(lm_late)[2]
        B_Dose_po <- exp(coef(lm_late)[1])

        #residuals for early data
        po_early$logresid <- log(po_early$ValueDose) -
          predict(lm_late, newdata = po_early)
        #linear regression of residuals gives alpha, A/Dose
        lm_early_resid <- lm(logresid ~ Time,
                       data = po_early)
        alpha_po <- -coef(lm_resid_early)[2]
        A_Dose_po <- exp(coef(lm_resid_early)[1])

        #residuals for absorption phase
        po_abs$logresid <- log(po_abs$ValueDose) -
          predict(lm_late, newdata = po_abs) -
          predict(lm_resid_early, newdata = po_abs)
        #linear regression of residuals gives kgutabs
        lm_resid_abs <- lm(logresid ~ Time,
                       data = po_abs)
        kgutabs_po <- -coef(lm_resid_abs)[2]

        #get k21, k12, kel from A, B, alpha, beta
        #source: https://www.boomer.org/c/p4/c19/c1903.php
        k21_po <- (A_Dose_po * beta_po + B_Dose_po*alpha_po) / (A_dose_po + B_Dose_po)
        kel_po <- (alpha_po * beta_po) / k21_po
        k12_po <- alpha_po + beta_po - k21_po - kel_po

        #Solve A_Dose for Fgutabs/V1
        #See cp_2comp() for equation for A that was solved for Fgutabs/V1
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
                                 msg = "Avg of method of residuals on PO and IV data")
          k12_avg <- mean(c(k12_iv,
                            k12_po),
                          na.rm = TRUE)
          par_DF <- assign_start(param_name = "k12",
                                 param_value = k12_avg,
                                 par_DF = par_DF,
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


    #sigma: Assign starting value as the median concentration
  TRYSIGMA <- stats::median(fitdata$Value, na.rm = TRUE)
  par_DF <- assign_start(param_name = grep(x = par_DF$param_name,
                                           pattern= "sigma",
                                           value = TRUE),
                         param_value = TRYSIGMA,
                         msg = "Median concentration",
                         par_DF = par_DF,
                         start_from = start_from_data)

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
  lower <- par_DF[param_name, "lower_bound"]
  upper <- par_DF[param_name, "upper_bound"]

  if(all(is.null(start_from) |
         is.na(start_from) |
         length(start_from)==0 |
         !nzchar(start_from))){
    start_from <- NULL
  }

  if(all(tolower(start_from) %in% "all")){
    start_from <- par_DF$param_name
  }

  #test bounds
  if((param_value >= lower) %in% TRUE &
     (param_value <= upper) %in% TRUE &
     param_name %in% start_from){     #if within bounds, update par_DF
    par_DF[param_name, c("start_value",
                         "start_value_msg")] <- list(param_value,
                                                     msg)
  } #if outside bounds, return par_DF unchanged

  return(par_DF)
}
