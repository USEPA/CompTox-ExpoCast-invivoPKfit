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
#' | param_name     | start_value | start_value_msg | | ---------------|
#' ----------- | --------------- | | kelim          | 0.25        | Default
#' | | Vdist          | 5.56        | Default         | | kgutabs        | 2.19
#' | Default         | | Fgutabs        | 0.5         | Default         | | V1 |
#' 5.56           | Default         | | k12            | 0.2         | Default |
#' | k21            | 0.5         | Default         | | Fgutabs_Vdist  |
#' 0.5/5.56    | Default         | | Fgutabs_V1     | 0.5/5.56    | Default | |
#' sigma          | 1         | Default         |
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
#' A starting value for `A` is estimated as the median dose-normalized
#' concentration in `fitdata`.
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
#' In case there is not enough data available in each phase to perform
#' regression, the default or `httk`-derived starting values are used.
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
#' Then, the "early" data are predicted using the "late" linear regression, and
#' the residuals are calculated. A second linear regression is performed on the
#' "early" residuals. The negative slope of this second linear regression
#' provides an estimate for `kgutabs`. Then, the intercept of the late
#' regression, `A`, can be used with the estimates for `kelim` and `kgutabs` to
#' calculate an estimate for `Fgutabs_Vdist`:
#'
#' -  `Fgutabs_Vdist <- A * (kgutabs - kelim) / kgutabs`
#'
#' In case there is not enough data available in each phase to perform
#' regression, the default or `httk`-derived starting values are used.
#'
#' #### Both PO and IV dosing data available
#'
#' The final starting value for `kelim` is taken as the average of the IV-based
#' estimate and the PO-based estimate.
#'
#' The estimates for `Vdist` (from the IV data) and `Fgutabs_Vdist` (from the PO
#' data) are multiplied to produce an estimate for `Fgutabs`.
#'
#' In case there is not enough data available in each phase to perform
#' regression, the default or `httk`-derived starting values are used.
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
#' "elbow point" of the data (using [find_elbow()]).
#'
#' First, a linear regression is performed on the "late" data. The negative
#' slope of the "late" linear regression is recorded as `beta`, and the
#' exponentiated intercept as `B`.
#'
#' Then the "early" data are predicted using the "late" linear regression, and
#' the residuals are calculated. A second linear regression is performed on the
#' "early" residuals. The negative slope of this second linear regression is
#' recorded as `alpha`, and the exponentiated intercept as `A`. Then, the
#' parameter starting values are calculated as follows (see
#' https://www.boomer.org/c/p4/c19/c1902.php):
#'
#' -`k21 <- (A * beta + B*alpha)/(A + B)`
#' -`kelim <- (alpha * beta)/k21`
#' -`k12 <- alpha + beta - k21 - kel`
#' -`V1 <- (alpha - k21)/(A * (alpha - beta))`
#'
#' In case there is not enough data available in each phase to perform
#' regression, the default or `httk`-derived starting values are used.
#'
#' ### PO dosing data available
#'
#' `kelim`, `k12`, `k21` `kgutabs`, and `Fgutabs_V1` (the ratio of `Fgutabs` to
#' `V1`) are estimated using the method of residuals on `log(Value/Dose)` vs.
#' `Time` from the PO data. Specifically, the data are divided into three parts:
#' "absorption" phase, "early" phase, and "late" phase. The time of peak
#' concentration is taken as the end of the "absorption" phase data, found using
#' [find_peak()]. Then, for all data after the time of peak concentration, the
#' "elbow" point is found (using [find_elbow()]) and used as the dividing line
#' between early and late-phase data.
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
#' Using `alpha`, `A`, `beta`, and `B`, estimates for `k21`, `kelim`, and `k12`
#' can be made:
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
#' In case there is not enough data available in each phase to perform
#' regression, the default or `httk`-derived starting values are used.
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
#' In case there is not enough data available in each phase to perform
#' regression, the default or `httk`-derived starting values are used.
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
#'   estimated from the data, so long as the `httk`-derived parameter values are
#'   within the lower and upper bounds defined by `par_DF$lower` and
#'   `par_DF$upper`. To suppress estimation of starting values from data, set
#'   `start_from_data` to NULL, NA, or a zero-length character vector. If a
#'   parameter name is in both `start_from_httk` and `start_from_data`, then
#'   data-derived starting values will override `httk` starting values for that
#'   parameter.
#' @param fit_log_conc Logical: Whether you are fitting log-transformed
#'   concentrations, or not. If TRUE, then `sigma` represents the log-scale
#'   standard devation, and a starting value for `sigma` will be determined
#'   using the residual error between log-scale observed and predicted
#'   concentrations. If FALSE (default), then `sigma` represents the
#'   natural-scale residual standard deviation, and a starting value for `sigma`
#'   will be calculated using the residual error between natural-scale predicted and observed
#'   concentrations.
#' @param suppress.messages Logical: Whether to suppress printing of messages.
#'
#' @return The data.frame `par_DF` with added variables `start_value`,
#'   containing each parameter's starting value, and `start_value_msg`,
#'   containing a brief message about how each starting value was calculated.
#' @export
#' @import data.table
#' @author Caroline Ring, John Wambaugh, Mitchell Teague

get_starts <- function(par_DF = NULL,
                       model,
                       fitdata,
                       pool_sigma = FALSE,
                       fit_conc_dose = TRUE,
                       starts_default = data.frame(
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
                         start_value = c(0.25, #kelim
                                         5.56, #Vdist
                                         2.19, #kgutabs
                                         0.5, #Fgutabs
                                         5.56, #V1
                                         0.2, #k12
                                         0.5, #k21
                                         0.5/5.56, #Fgutabs_Vdist
                                         0.51/5.56, #Fgutabs_V1
                                         1, #Rblood2plasma
                                         1), #sigma
                         start_value_msg = "Default"
                       ),
                       start_from_httk = "all",
                       start_from_data = "all",
                       fit_log_conc = FALSE,
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

    #the abovemerge won't get study-specific sigmas: handle them
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
    } #end for loop over study-specific sigmas

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

        #set "Fgutabs_Vdist" or "Fgutabs_V1" appropriately
        httk_params[["Fgutabs_Vdist"]] <- httk_params[["Fgutabs"]]/httk_params[["Vdist"]]
        httk_params[["Fgutabs_V1"]] <- httk_params[["Fgutabs"]]/httk_params[["V1"]]

      #replace the parameters in par_DF whose names match those in httk_params,
      #and are named in start_from_Httk
      #find matching names
      replace_names <- intersect(names(httk_params),
                                       intersect(par_DF$param_name,
                                                 start_from_httk))

      #if optimize_params is FALSE for Rblood2plasma, do not use its httk value
      #keep the default of 1
      if("Rblood2plasma" %in% par_DF$param_name){
      if(par_DF[par_DF$param_name %in% "Rblood2plasma",
                "optimize_param"] %in% FALSE){
        replace_names <- setdiff(replace_names,
                                 "Rblood2plasma")
      }
      }

      #update par_DF for any of these parameters
      for (this_param in replace_names){
        if(this_param %in% par_DF$param_name){
          par_DF <- assign_start(param_name = this_param,
                                 param_value = httk_params[[this_param]],
                                 msg = paste0("httk::parameterize_1comp(",
                                              "dtxsid = '", this.dtxsid, "', ",
                                              "default.to.human = TRUE, ",
                                              "suppress.messages = TRUE, ",
                                              "species = '", httk_species,
                                              "')"),
                                 start_from = start_from_httk,
                                 par_DF <- par_DF)
        }
      }


    } #end if(this.dtxsid %in% httk::get_cheminfo(info="dtxsid"))

    #Try roughly estimating parameters from data

    #Drop control points & non-detects
    tmpdat <- subset(fitdata,
                     Dose > 0)

    tmpdat <- subset(tmpdat,
                     is.finite(Value))

    #log transform Value/Dose
    tmpdat$logValueDose <- log(tmpdat$Value_Dose)

    #Split into IV and oral datasets
    iv_data <- tmpdat[tmpdat$Route %in% "iv", ] #will be empty if no IV data
    po_data <- tmpdat[tmpdat$Route %in% "po", ] #will be empty if no PO data

    has_iv <- any(tmpdat$Route %in% "iv")
    has_po <- any(tmpdat$Route %in% "po")

    has_blood <- any(tmpdat$Media %in% "blood")
    has_plasma <- any(tmpdat$Media %in% "plasma")

    #flat model
    if(model %in% "flat"){
      Rblood2plasma_iv <- NA_real_
      Rblood2plasma_po <- NA_real_

      if(has_iv %in% TRUE){
        if(has_plasma %in% TRUE){
          Vdist_log <- -1 * mean(iv_data[iv_data$Media %in% "plasma",
                                          "logValueDose"],
                                  na.rm = TRUE)
          Vdist <- exp(Vdist_log)
          if(has_blood %in% TRUE){ #if both blood and plasma IV data, estimate Rblood2plasma
            ymean <- mean(iv_data[iv_data$Media %in% "blood",
                                                  "logValueDose"],
                                          na.rm = TRUE)
            Rblood2plasma_iv_log <- ymean + Vdist_log
            Rblood2plasma_iv <- exp(Rblood2plasma_iv_log)
          }
        }else{ #if blood data only
          #just estimate a blood Vdist
          Vdist_log <- mean(iv_data[iv_data$Media %in% "blood",
                                             "logValueDose"],
                                      na.rm = TRUE)
          Vdist <- exp(Vdist_log)
        }

        par_DF <- assign_start(param_name = "Vdist",
                               param_value = Vdist,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Based on average Value/Dose from IV data")
      }

      if(has_po %in% TRUE){
        if(has_plasma %in% TRUE){
        Fgutabs_Vdist_log <- mean(po_data[po_data$Media %in% "plasma",
                                                      "logValueDose"],
                                              na.rm = TRUE)
        Fgutabs_Vdist <- exp(Fgutabs_Vdist_log)
        if(has_blood %in% TRUE){ #if both blood and plasma IV data, estimate Rblood2plasma
          ymean <- mean(po_data[po_data$Media %in% "blood",
                                "logValueDose"],
                        na.rm = TRUE)
          Rblood2plasma_po_log <- ymean - Fgutabs_Vdist_log
          Rblood2plasma_po <- exp(Rblood2plasma_po_log)
        }
        }else{
          #if blood-only data, estimate Fgutabs_Vdist for blood
          Fgutabs_Vdist_log <- mean(po_data[po_data$Media %in% "blood",
                                                   "logValueDose"],
                                           na.rm = TRUE)
          Fgutabs_Vdist <- exp(Fgutabs_Vdist_log)
        }

        #if both oral and IV data, then use Vdist from IV to calculate Fgutabs
        if(has_iv %in% TRUE){
          Fgutabs <- Fgutabs_Vdist * Vdist
          par_DF <- assign_start(param_name = "Fgutabs",
                                 param_value = Fgutabs,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Based on average Value/Dose from PO and IV data")
        }

        par_DF <- assign_start(param_name = "Fgutabs_Vdist",
                               param_value = Fgutabs_Vdist,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Based on average Value/Dose from PO data")
      }

      if(has_blood %in% TRUE & has_plasma %in% TRUE){
      Rblood2plasma <- mean(c(Rblood2plasma_iv,
                              Rblood2plasma_po),
                            na.rm = TRUE)

      par_DF <- assign_start(param_name = "Rblood2plasma",
                             param_value = Rblood2plasma,
                             par_DF = par_DF,
                             start_from = start_from_data,
                             msg = "Based on comparing average Value/Dose from blood and plasma data")
      }

    }else{ #for non-flat models



    #####################
    # 1-compartment model
    #####################
    if(model %in% "1compartment"){
      if(has_iv %in% TRUE){
        #--------------------------------------------------
        #get kelim and Vdist from linear regression
        #see https://www.boomer.org/c/p4/c04/c0406.php
        #--------------------------------------------------
        lm_iv <- do_linreg(x = iv_data$Time,
                           y = iv_data$logValueDose,
                           intercept_exp = TRUE,
                           slope_neg = TRUE)
        Vdist_iv <- 1/lm_iv$intercept
        kelim_iv <- lm_iv$slope

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

        po_peak <- get_peak(x = po_data$Time,
                            y = po_data$logValueDose)

        po_early <- subset(po_data,
                           Time <= po_peak$x)

        po_late <- subset(po_data,
                          Time >= po_peak$x)

        A_Dose_po <- NA_real_
        kelim_po <- NA_real_
        kgutabs_po <- NA_real_
        Fgutabs_Vdist_po <- NA_real_
        kelim_po_msg <- NA_character_
        kgutabs_po_msg <- NA_character_
        Fgutabs_Vdist_po_msg <- NA_character_

        #if we have an IV estimate of kelim, use it
        have_kelim_iv <- par_DF[par_DF$param_name %in% "kelim",
                                "start_value_msg"] %in%
          "Linear regression of IV data"
        if(have_kelim_iv %in% TRUE){
          kelim_po <- kelim_iv
          kelim_po_msg <- "Linear regression on IV data"
        }else{ #if we do not have an IV estimate of kelim, get one
        if(length(unique(po_late$Time))>2){ #if we have any elimination-phase data
          #try to estimate kelim via linear regression
          #keep the latter half of po_late, or enough to do regression
          po_late_timept <- unique(po_late$Time)
          po_late_timept <- sort(po_late_timept, decreasing =TRUE)
          #midpoint of time points after peak
          late_midpt <- mean(c(po_peak$x, max(po_late_timept)))
          #keep timepoints after midpoint, or at least 3 timepoints
          late_keep <- pmax(sum(po_late_timept >= late_midpt),
                            3)
          po_late <- po_late[po_late$Time <= max(po_late_timept[1:late_keep]), ]

        lm_po_late <- do_linreg(x = po_late$Time,
                                y = po_late$logValueDose,
                                intercept_exp = TRUE,
                                slope_neg = TRUE)
        A_Dose_po <- lm_po_late$intercept
        kelim_po <- lm_po_late$slope
        kelim_po_msg <- "Linear regression on late-phase PO data"

        if(!is.finite(kelim_po)){ #if regression didn't produce a valid estimate
         #make a gross approximation
            #assume tmax = 5 *elimination half-life
          thalf_elim <- max(po_data$Time)/5
          kelim_po <- log(2)/thalf_elim
          kelim_po_msg <- "Assume PO data tmax = 5 * elimination half-life"
          }

        #Assign kelim starting value
        par_DF <- assign_start(param_name = "kelim",
                               param_value = kelim_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = kelim_po_msg)
        }else{
          #if we only have absorption phase, not elimination phase
          #then make a gross guess:
          #assume tpeak is 2 * absorption halflife
          thalf_abs <- po_peak$x/2
          kgutabs_po <- log(2)/thalf_abs
          kgutabs_po_msg <- paste("No elimination-phase data:",
          "fewer than two detects after tpeak.",
          "Assume tpeak = 2*absorption halflife.")
          par_DF <- assign_start(param_name = "kgutabs",
                                 param_value = kgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = kgutabs_po_msg)
          #another gross guess: kelim is same as kgutabs
          #(unless we have IV estimate of kelim)
          kelim_po <- ifelse(is.finite(kelim_iv),
                             kelim_iv,
                             kgutabs_po/10)
          kelim_po_msg <- ifelse(is.finite(kelim_iv),
                                 "Linear regression of IV data",
                                 paste("No elimination-phase data:",
                                       "fewer than two detects after tpeak.",
                                 "Assume tpeak = 2*absorption halflife",
                                 "and kelim = kgutabs/10"))
          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kelim_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = kelim_po_msg)
          }
        } #end if(have_kelim_iv %in% TRUE)/else



#Continue with approximations of kgutabs

#First, try Wagner-Nelson method
        #see https://www.boomer.org/c/p4/c09/c0903.php
        #get AUC/Dose

        #sort data by dose, then by increasing time
        po_data <- as.data.table(po_data)

        setorder(po_data,
                 Dose,
                 Time)

        #get average Value for every dose & time point
        #weight by number of subjects
        po_avg <- po_data[, .(
          Value_avg =exp(
            weighted.mean(
              log(Value),
              N_Subjects)
            )
          ),
                by = .(Dose, Time)]

        setorder(po_avg, Dose, Time)

        #calculate cumulative AUCs using trapezoidal rule
        #assume everything starts at (0,0)
        po_avg[, AUC := pracma::cumtrapz(x = c(0, Time),
                                          y = c(0,Value_avg)
                                         )[2:(.N+1), 1],
               by = Dose]
        #Calculate A/V
        po_avg[, A_V := Value_avg + kelim_po * AUC]

        #To calculate Amax/V: First calculate AUC_infinity
        po_avg[, AUC_inf := AUC[.N] + Value_avg[.N] / kelim_po,
               by = Dose]

        #Calculate Amax/V
        po_avg[, Amax_V := kelim_po * AUC_inf]

        #Calculate Amax/V - A/V
        po_avg[, Amax_minus_A := Amax_V - A_V]

        #keep only Amax - A values that can be log-transformed (positive)
        po_avg <- po_avg[Amax_minus_A > 0, ]

        #Normalize by Dose
        po_avg[, Amax_minus_A_Dose := Amax_minus_A/Dose]

        setorder(po_avg, Time)

        #take the early points:
        #those before po_peak$x,
        #and enough extras to allow linear regression, if necessary
        n_keep <- pmax(po_avg[Time <= po_peak$x, .N],
                      3)
        po_avg_early <- po_avg[seq(1, n_keep), ]

        #method of Wagner-Nelson:
        #linear regression of log((Amax - A)/V) vs. time
        #should have slope -kgutabs and intercept F*Dose/V
        #When dose-normalized, the intercept will be F/V
        wagnel <- do_linreg(x = po_avg_early$Time,
                  y = log(po_avg_early$Amax_minus_A_Dose),
                  slope_neg = TRUE,
                  intercept_exp = TRUE)

        kgutabs_po <- wagnel$slope
        kgutabs_po_msg <- "Wagner-Nelson method on early PO data"

        Fgutabs_Vdist_po <- wagnel$intercept
        Fgutabs_Vdist_po_msg <- "Wagner-Nelson method on early PO data"

if(is.finite(kgutabs_po) &
   is.finite(Fgutabs_Vdist_po)){
        par_DF <- assign_start(param_name = "kgutabs",
                               param_value = kgutabs_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = kgutabs_po_msg)
        par_DF <- assign_start(param_name = "Fgutabs_Vdist",
                               param_value = Fgutabs_Vdist_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = Fgutabs_Vdist_po_msg)

}else{ #if Wagner-Nelson failed, try method of residuals to get kgutabs
  #for this, we need to have >1 points in early-phase data
  #if we did not get an intercept from linear regression of late-phase data,
  # back-calculate it from the kelim estimate
  #and the peak time/ & peak conc
  if(!is.finite(A_Dose_po)){
    A_Dose_po <- exp(po_peak$y +
                       kelim_po * po_peak$x)
  }

  #get residuals
  #predicted values:
  pred_late <- exp(log(A_Dose_po) - kelim_po * po_early$Time)
  #residuals: late predictions - early observations
  #because late predictions should be > early observations
  resid_early <- pred_late - po_early$Value_Dose
  #log-transform residuals
  suppressWarnings(po_early$logresid <- log(resid_early))
  #drop any NA log residuals
  po_early <- subset(po_early, is.finite(logresid))

  lm_resid_early <- do_linreg(x = po_early$Time,
                              y = po_early$logresid,
                              intercept_exp = TRUE,
                              slope_neg = TRUE)
  kgutabs_po <- lm_resid_early$slope

  #update par_DF
  #kgutabs
  kgutabs_po_msg <- "Method of residuals on PO data"
  par_DF <- assign_start(param_name = "kgutabs",
                         param_value = kgutabs_po,
                         par_DF = par_DF,
                         start_from = start_from_data,
                         msg = kgutabs_po_msg)

  #and get Fgutabs_Vdist_po from intercept
  Fgutabs_Vdist_po <- A_Dose_po * (kgutabs_po - kelim_po)/kgutabs_po
  Fgutabs_Vdist_po_msg <- "Calc from intercept of regression on late-phase data"
  par_DF <- assign_start(param_name = "Fgutabs_Vdist",
                         param_value = Fgutabs_Vdist_po,
                         par_DF = par_DF,
                         start_from = start_from_data,
                         msg = Fgutabs_Vdist_po_msg)
}

#if both Wagner-Nelson and method of residuals failed to get kgutabs
        if(!is.finite(kgutabs_po)){
        #Make a gross assumption that tpeak = 2 *
        #absorption half-life -- see https://www.boomer.org/c/p4/c09/c0904.php
          thalf_abs <- po_peak$x/2
          kgutabs_po <- log(2)/thalf_abs
          #update par_DF
          #kgutabs
          kgutabs_po_msg <- paste("Wagner-Nelson and method of residuals",
          "on PO data both failed.",
          "Assume tpeak = 2 * absorp half-life")
          par_DF <- assign_start(param_name = "kgutabs",
                                 param_value = kgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = kgutabs_po_msg)

          if(!is.finite(A_Dose_po)){
            A_Dose_po <- exp(po_peak$y +
                               kelim_po * po_peak$x)
          }

          #and get Fgutabs_Vdist_po from intercept
          Fgutabs_Vdist_po <- A_Dose_po * (kgutabs_po - kelim_po)/kgutabs_po
          Fgutabs_Vdist_po_msg <- "Calc from intercept of regression on late-phase data using estimates of kgutabs and kelim"
          par_DF <- assign_start(param_name = "Fgutabs_Vdist",
                                 param_value = Fgutabs_Vdist_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = Fgutabs_Vdist_po_msg)

        }

        #Fgutabs_Vdist

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
                                 msg = paste0(
                                   "Method of residuals on PO data to get Fgutabs/Vdist (",
                                   Fgutabs_Vdist_po, ")",
                                   " combined with estimate of Vdist")
          )


        }
      } #end if(has_po %in% TRUE)
    } #end if(model %in% "1compartment")

    #####################
    # 2-compartment model
    #####################
    if(model %in% "2compartment"){
      if(has_iv %in% TRUE){
        #--------------------------------------------------
        #get alpha, beta, A, and B from method of residuals
        #see https://www.boomer.org/c/p4/c19/c1903.php
        #--------------------------------------------------

        #split IV data into early and late parts
        #find the dividing line between early and late as the "elbow" point
        elbow <- get_elbow(x = iv_data$Time,
                           y = iv_data$logValueDose)

        #split IV data into early and late phases at elbow_time
        iv_early <- subset(iv_data,
                           Time <= elbow$x)
        iv_late <- subset(iv_data,
                          Time >= elbow$x)

        #--------------------------------------------------
        # Regression on late-phase data
        #-------------------------------------------------

        lm_late <- do_linreg(x = iv_late$Time,
                             y = iv_late$logValueDose,
                             intercept_exp = TRUE,
                             slope_neg = TRUE)
        B_Dose_iv <- lm_late$intercept
        beta_iv <- lm_late$slope

        #if regression failed for beta, assume max time is 2 * terminal half-life
        if(!is.finite(beta_iv) | beta_iv <= 0){
          #use max time for all IV data, whether detect or no
          iv_all <- subset(fitdata, Route %in% "iv")
          thalf_beta <- max(iv_all$Time)/2
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
        resid_early <- iv_early$Value_Dose - pred_late
        suppressWarnings(iv_early$logresid <-  log(resid_early))
        #drop any NA residuals
        iv_early <- subset(iv_early, is.finite(logresid))

        lm_resid_early <- do_linreg(x = iv_early$Time,
                                     y = iv_early$logresid,
                                    intercept_exp = TRUE,
                                    slope_neg = TRUE)

        A_Dose_iv <- lm_resid_early$intercept
        alpha_iv <- lm_resid_early$slope

        #if method of residuals failed for alpha,
        #then try assuming that elbow point = 2 * early-phase half life
        if(!is.finite(alpha_iv) | alpha_iv <= 0){
          thalf_alpha_iv <- elbow$x/2
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
      }else{ #if !(has_iv %in% TRUE)
        k21_iv <- NA_real_
        kel_iv <- NA_real_
        k12_iv <- NA_real_
        V1_iv <- NA_real_
      }

      #----------------------------------------------------------------
      # Oral data -- 2 compartment model
      # See https://www.boomer.org/c/p4/c19/c1907.php
      # Use method of residuals *twice*
      #----------------------------------------------------------------

      if(has_po %in% TRUE){
        #Split the data into three parts:
        #Absorption phase, early phase, and late phase.

        #Absorption phase ends when peak concentration is achieved:
        po_peak <- get_peak(x = po_data$Time,
                            y = po_data$logvalueDose)

        #split data into absorption and non-absorption phases
        po_abs <- subset(po_data, Time <= po_peak$x)
        po_nonabs <- subset(po_data, Time >= po_peak$y)

        #further split the non-absorption phase into early and late parts
        #find the dividing line between early and late as the "elbow" point
        elbow <- get_elbow(x = po_nonabs$Time,
                           y = po_nonabs$logValueDose)

        po_early <- subset(po_nonabs, Time <= elbow$x)
        po_late <- subset(po_nonabs, Time >= elbow$x)

        lm_late <- do_linreg(x = po_late$Time,
                             y = po_late$logValueDose,
                             intercept_exp = TRUE,
                             slope_neg = TRUE)

        B_Dose_po <- lm_late$intercept
        beta_po <- lm_late$slope

        B_Dose_po_fail <- FALSE

        #if regression fails for beta_po,
        #assume last time point is 2 times terminal half-life
        if(!is.finite(beta_po) |
           beta_po < 0){
          #use max time for all po data, whether detect or nondetect
          po_all <- subset(fitdata, Route == "po")
          thalf_beta <- max(po_all$Time)/2
          beta_po <- log(2)/thalf_beta
          #intercept: extrapolate back to time = 0 from elbow point
          B_Dose_po <- exp(elbow$y + beta_po * elbow$x)
          B_Dose_po_fail <- TRUE
        }

        #residuals for early data
        #obs should be greater than predicted
        #calc residuals during early phase observed - (predicted from late-phase
        #regression) (because late-phase predictions should be less than
        #early-phase observations)
        pred_late_po <- exp(log(B_Dose_po) + -beta_po * po_early$Time)
        resid_early_po <- po_early$Value_Dose - pred_late_po
        suppressWarnings(po_early$logresid <- log(resid_early_po))
      #drop any NA log residuals
        po_early <- subset(po_early,
                           is.finite(logresid))

        lm_resid_early <- do_linreg(x = po_early$Time,
                                    y = po_early$logresid,
                                    intercept_exp = TRUE,
                                    slope_neg = TRUE)
        A_Dose_po <- lm_resid_early$intercept
        alpha_po <- lm_resid_early$slope
        A_Dose_po_fail <- FALSE

        #if method of residuals failed for alpha,
        #then try assuming that elbow point = 2 * early-phase half life
        if(!is.finite(alpha_po) | alpha_po <= 0){
          thalf_alpha <- elbow$x/2
          alpha_po <- log(2)/thalf_alpha
          #extrapolate from elbow point back to time = 0
            A_Dose_po <- exp(elbow$y + alpha_po * elbow$x)
            A_Dose_po_fail <- TRUE
        }

        #residuals for absorption phase
        #predictions from late-phase regression:
        late_pred <- exp(log(B_Dose_po) + -beta_po * po_abs$Time)
        #predictions from early-phase residuals regression
        early_pred <- exp(log(A_Dose_po) + -alpha_po * po_abs$Time)

        #predictions should be > observations now, so reverse sign to calc resids
        po_abs$resid <- (late_pred + early_pred) - po_abs$Value_Dose
        suppressWarnings(po_abs$logresid <- log(po_abs$resid))

        #drop NA residuals
        po_abs <- subset(po_abs,
                         is.finite(logresid))

        lm_resid_abs <- do_linreg(x = po_abs$Time,
                                  y = po_abs$logresid,
                                  intercept_exp = TRUE,
                                  slope_neg = TRUE)
        kgutabs_po <- lm_resid_abs$slope

        par_DF <- assign_start(param_name = "kgutabs",
                               param_value = kgutabs_po,
                               par_DF = par_DF,
                               start_from = start_from_data,
                               msg = "Method of residuals on PO data")

        #if method of residuals failed, then try just assuming that tpeak = 5 *
        #absorption half-life -- see https://www.boomer.org/c/p4/c09/c0904.php
        if(!is.finite(kgutabs_po) | kgutabs_po <= 0){
          thalf_abs <- po_peak$x/2
          kgutabs_po <- log(2)/thalf_abs
          par_DF <- assign_start(param_name = "kgutabs",
                                 param_value = kgutabs_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of inspection on PO data, assume tpeak = 5 * absorption half-life")
        }

        #get k21, k12, kel from A, B, alpha, beta. Source: derived from formulas
        #in spreadsheet at https://www.boomer.org/c/p4/c19/c1907.php. NB: PO A
        #and B are not the same as IV A and B, so the formula for PO k21 is not the
        #same as the formula for IV k21.
        k21_po <- (A_Dose_po * (kgutabs_po - alpha_po) * beta_po +
                     B_Dose_po * (kgutabs_po - beta_po) * alpha_po) /
          (A_Dose_po* (kgutabs_po - alpha_po) +
             B_Dose_po * (kgutabs_po - beta_po))

        kel_po <- (alpha_po * beta_po) / k21_po
        k12_po <- alpha_po + beta_po - k21_po - kel_po

        #Solve A_Dose for Fgutabs/V1
        #See cp_2comp()
        Fgutabs_V1_po_A <- A_Dose_po *
          ( (kgutabs_po - alpha_po) * (beta_po - alpha_po))/
          ( kgutabs_po * (k21_po - alpha_po) )
        #Solve B_Dose for Fgutabs/V1
        #See cp_2comp()
        Fgutabs_V1_po_B <- B_Dose_po *
          ( (kgutabs_po - beta_po) * (alpha_po - beta_po)) /
          ( kgutabs_po * (k21_po - beta_po))

        if(A_Dose_po_fail %in% TRUE &
           B_Dose_po_fail %in% FALSE){
          #if late regression was OK but early regression failed,
          #use estimate from late regression intercept
          Fgutabs_V1_po <- Fgutabs_V1_po_B
        }else{
          #just take the average from early and late regression intercepts
          Fgutabs_V1_po <- mean(c(Fgutabs_V1_po_A,
                                 Fgutabs_V1_po_B),
                                 na.rm = TRUE)
        }

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
        }

          #k21, kelim, k12:
          #update par_DF using PO estimates only if no IV data,
        #or if IV fitting failed to produce an acceptable estimate
         if(par_DF[par_DF$param_name %in% "k21",
                   "start_value_msg"] != "Method of residuals on IV data"){
          par_DF <- assign_start(param_name = "k21",
                                 param_value = k21_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data"
                                 )
         }

          if(par_DF[par_DF$param_name %in% "kelim",
                    "start_value_msg"] != "Method of residuals on IV data"){
          par_DF <- assign_start(param_name = "kelim",
                                 param_value = kel_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data"
                                 )
          }

          if(par_DF[par_DF$param_name %in% "k12",
                    "start_value_msg"] != "Method of residuals on IV data"){
          par_DF <- assign_start(param_name = "k12",
                                 param_value = k12_po,
                                 par_DF = par_DF,
                                 start_from = start_from_data,
                                 msg = "Method of residuals on PO data")
          }
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
                       iv.dose = fitdata$Route %in% "iv",
                       medium = fitdata$Media))

  if(fit_conc_dose %in% TRUE){
    if(fit_log_conc %in% FALSE){
    resid <- pred/fitdata$Dose - fitdata$Value_Dose
    }else{
      resid <- log(pred/fitdata$Dose) - log(fitdata$Value_Dose)
    }
  }else{
    if(fit_log_conc %in% FALSE){
    resid <- pred - fitdata$Value
    }else{
      resid <- log(pred) - log(fitdata$Value)
    }
  }


  resid[!is.finite(resid)] <- NA_real_

  for(this_sigma in grep(x = par_DF$param_name,
                         pattern= "sigma",
                         value = TRUE)){
    if(this_sigma == "sigma"){ #if only one sigma (one study, or pooled)
      tmp_sigma <- sd(resid, na.rm = TRUE)
    }else{
      #if multiple sigmas for multiple studys
      #get the study
      studyid <- gsub(x = this_sigma,
                    pattern = "sigma_study_",
                    replacement = "")
      resid_study <- resid[fitdata$Study %in% studyid]
      tmp_sigma <- sd(resid_study, na.rm = TRUE)
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
     (param_name %in% start_from) %in% TRUE){
    if((param_value >= lower) %in% TRUE &
     (param_value <= upper) %in% TRUE){     #if within bounds and finite, update par_DF
    par_DF[par_DF$param_name %in% param_name, c("start_value",
                         "start_value_msg")] <- list(param_value,
                                                     msg)
     }else if((param_value >= lower) %in% FALSE){
    #use the lower bound
       lb <- par_DF[par_DF$param_name %in% param_name, "lower_bound"]
       par_DF[par_DF$param_name %in% param_name,
              c("start_value",
                                                   "start_value_msg")] <- list(lb,
                                                                               paste("The following estimated a value below the lower bound (",
                                                                                     signif(param_value,3),
                                                                               "); the lower bound was substituted. Original message:",
                                                                                     msg))
     }else if((param_value <= upper) %in% FALSE){
    #use the upper bound
       ub <- par_DF[par_DF$param_name %in% param_name, "upper_bound"]
       par_DF[par_DF$param_name %in% param_name,
              c("start_value",
                "start_value_msg")] <- list(ub,
                                            paste("The following estimated a value above the upper bound (",
                                                  signif(param_value,3),
                                                  "); the upper bound was substituted. Original message:",
                                                  msg))
  }
    #if not finite, return par_DF unchanged
  }
  }

  return(par_DF)
}
