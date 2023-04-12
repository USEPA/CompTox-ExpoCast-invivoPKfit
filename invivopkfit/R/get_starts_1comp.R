get_starts_1comp <- function(data,
                             par_DF){

  #Drop control points & non-detects
  tmpdat <- subset(data,
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
    #get kelim and Fgutabs_Vdist from linear regression
    #see https://www.boomer.org/c/p4/c09/c0901.php
    #--------------------------------------------------
    #split data into absorption and elimination phases,
    #using the time of peak concentration as the dividing line

    po_peak <- get_peak(x = po_data$Time,
                        y = po_data$logValueDose)

    po_abs <- subset(po_data,
                       Time <= po_peak$x)

    po_elim <- subset(po_data,
                      Time >= po_peak$x)
#initialize params to be estimated via regression
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
      if(length(unique(po_elim$Time))>2){ #if we have any elimination-phase data
        #try to estimate kelim via linear regression
        #keep the latter half of po_elim, or enough to do regression
        #get time points during elimination phase, sorted in decreasing order
        po_elim_timept <- unique(po_elim$Time)
        po_elim_timept <- sort(po_elim_timept, decreasing =TRUE)
        #find midpoint of time points after peak
        late_midpt <- mean(c(po_peak$x, max(po_elim_timept)))
        #keep timepoints after midpoint, or at least 3 timepoints, whichever is greater
        late_keep <- pmax(sum(po_elim_timept >= late_midpt),
                          3)
        po_elim <- po_elim[po_elim$Time <= max(po_elim_timept[1:late_keep]), ]
        #do linear regression
        lm_po_elim <- do_linreg(x = po_elim$Time,
                                y = po_elim$logValueDose,
                                intercept_exp = TRUE,
                                slope_neg = TRUE)
        A_Dose_po <- lm_po_elim$intercept
        kelim_po <- lm_po_elim$slope
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
    po_avg_abs <- po_avg[seq(1, n_keep), ]

    #method of Wagner-Nelson:
    #linear regression of log((Amax - A)/V) vs. time
    #should have slope -kgutabs and intercept F*Dose/V
    #When dose-normalized, the intercept will be F/V
    wagnel <- do_linreg(x = po_avg_abs$Time,
                        y = log(po_avg_abs$Amax_minus_A_Dose),
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
      pred_elim <- exp(log(A_Dose_po) - kelim_po * po_abs$Time)
      #residuals: late predictions - early observations
      #because late predictions should be > early observations
      resid_abs <- pred_elim - po_abs$Value_Dose
      #log-transform residuals
      suppressWarnings(po_abs$logresid <- log(resid_abs))
      #drop any NA log residuals
      po_abs <- subset(po_abs, is.finite(logresid))

      lm_resid_abs <- do_linreg(x = po_abs$Time,
                                  y = po_abs$logresid,
                                  intercept_exp = TRUE,
                                  slope_neg = TRUE)
      kgutabs_po <- lm_resid_abs$slope

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

  return(par_DF)
}
