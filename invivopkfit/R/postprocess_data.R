#' Postprocess fitted PK models
#'
#' Derive quantities such as half-life, tpeak, Cpeak, Css from fitted PK models
#'
#' @param PK_fit A data.table of fitted PK parameters, as produced by `fit_all()` or `analyze_subset()`
#' @param model The name of the relevant model: "1compartment", "2compartment", or "flat"
#' @return A data.table of derived quantities relevant to the model, with one row for each set of model parameters
#' @export postprocess_data
#' @author John Wambaugh, Caroline Ring
postprocess_data <- function(PK_fit,
                             model){

#Do post fit calculations

  if(model %in% "flat"){
    #reshape wide (one column for each parameter)
    PK_flat <- dcast(PK_fit[model %in% "flat" &
                               !grepl(x = param_name,
                                      pattern = "sigma")],
                      Analysis_Type + DTXSID + Species +
                        model + References.Analyzed +
                        Routes + Media +
                       time_units_fitted +
                        AIC ~ param_name,
                      value.var = "Fitted mean")
    #also get the SD for each fitted param
    PK_flat_sd <- dcast(PK_fit[model %in% "flat" &
                                  !grepl(x = param_name,
                                         pattern = "sigma")],
                         Analysis_Type + DTXSID + Species +
                           model + References.Analyzed +
                           Routes + Media +
                          time_units_fitted +
                           AIC ~ param_name,
                         value.var = "Fitted std dev")
    #append "_sd" for these params
    setnames(PK_flat_sd,
             PK_fit[model %in% "flat" &
                      !grepl(x = param_name,
                             pattern = "sigma"), unique(param_name)],
             PK_fit[model %in% "flat" &
                      !grepl(x = param_name,
                             pattern = "sigma"), paste0(unique(param_name),
                                                        "_sd")])
    #merge
    PK_flat <- merge(PK_flat,
                      PK_flat_sd,
                      by = intersect(names(PK_flat),
                                     names(PK_flat_sd))
    )

    PK_out <- copy(PK_flat)
  }else if(model %in% "1compartment"){
  #for 1-compartment model:

  #Calculate: total clearance; half-life; tpeak; Cpeak
  #and Css for a dose of 1 mg/kg/da

    #reshape wide (one column for each parameter)
  PK_1comp <- dcast(PK_fit[model %in% "1compartment" &
                             !grepl(x = param_name,
                                    pattern = "sigma")],
                    Analysis_Type + DTXSID + Species +
                      model + References.Analyzed +
                      time_units_fitted +
                      Routes + Media +
                      AIC ~ param_name,
                    value.var = "Fitted mean")
  #also get the SD for each fitted param
  PK_1comp_sd <- dcast(PK_fit[model %in% "1compartment" &
                                !grepl(x = param_name,
                                       pattern = "sigma")],
                    Analysis_Type + DTXSID + Species +
                      model + References.Analyzed +
                      Routes + Media +
                      time_units_fitted +
                      AIC ~ param_name,
                    value.var = "Fitted std dev")
  #append "_sd" for these params
  setnames(PK_1comp_sd,
           PK_fit[model %in% "1compartment" &
                    !grepl(x = param_name,
                           pattern = "sigma"), unique(param_name)],
           PK_fit[model %in% "1compartment" &
                    !grepl(x = param_name,
                           pattern = "sigma"), paste0(unique(param_name),
           "_sd")])
  #merge
 PK_1comp <- merge(PK_1comp,
                   PK_1comp_sd,
                   by = intersect(names(PK_1comp),
                                  names(PK_1comp_sd))
                   )

 #convert time constants back to hours if necessary
 time_units <- unique(PK_1comp$time_units_fitted)
 if(!(time_units %in% "hours")){
   time_const <- intersect(names(PK_1comp),
                           c("kelim",
                             "kgutabs",
                             "kelim_sd",
                             "kgutabs_sd")
   )
   if(time_units == "days"){
     convfun <- function(x) x*24
   }else if(time_units == "weeks"){
     convfun <- function(x) x*(24*7)
   }else if(time_units == "months"){
     convfun <- function(x) x*(24*30)
   }else if(time_units == "years"){
     convfun <- function(x) x*(24*365)
   }else if(time_units == "minutes"){
     convfun <- function(x) x/60
   }
   PK_1comp[, (time_const) := lapply(.SD,
                                     convfun),
            .SDcols = time_const]
 }

 PK_1comp[, time_units_reported := "hours"]

  PK_1comp[is.na(Fgutabs_Vdist),
           Fgutabs_Vdist := Fgutabs/Vdist]

  #CLtot:
  #If kelim and Vdist are available, CLtot = kelim * Vdist
  PK_1comp[, CLtot := kelim * Vdist]

  #if kelim and Fgutabs/Vdist are available, can only get CLtot/Fgutabs
  #kelim * (Vdist/Fgutabs) = CLtot/Fgutabs
  PK_1comp[, CLtot_Fgutabs := kelim / Fgutabs_Vdist]

  #Css_oral
  PK_1comp[, Css_oral_1mgkg := Fgutabs_Vdist * (1/(24*kelim))]
  #if only kelim and Vdist available -- get Css for IV infusion
  PK_1comp[, Css_iv_1mgkg := 1/(24 * CLtot)]

  #half-life, tpeak, Cpeak:

  #see https://www.boomer.org/c/p4/c08/c0803.php

  #half-life
  PK_1comp[, halflife := log(2) / kelim]

  #tpeak -- only available if kgutabs was fitted
  PK_1comp[, tpeak_oral := log(kgutabs / kelim) / (kgutabs - kelim)]

  #Cpeak for oral dose of 1 mg/kg
  PK_1comp[!is.na(kgutabs), Cpeak_oral_1mgkg := cp_1comp(params = list("kelim" = kelim,
                                                  "Fgutabs_Vdist" = Fgutabs_Vdist,
                                                  "kgutabs" = kgutabs),
                                    time = tpeak_oral,
                                    dose = 1,
                                    iv.dose = FALSE),
           by = .(Analysis_Type,
           DTXSID,
           Species,
             model,
           References.Analyzed,
             AIC)]

  PK_out <- copy(PK_1comp)
}else if(model %in% "2compartment"){

  #2-compartment model
  #reshape to wide format
  PK_2comp <- dcast(PK_fit[model %in% "2compartment" &
                             !grepl(x = param_name,
                                    pattern = "sigma")],
                    Analysis_Type + DTXSID + Species +
                      model + References.Analyzed +
                      Routes + Media +
                      time_units_fitted +
                      AIC ~ param_name,
                    value.var = "Fitted mean")

  #also get the SD for each fitted param
  PK_2comp_sd <- dcast(PK_fit[model %in% "2compartment" &
                                !grepl(x = param_name,
                                       pattern = "sigma")],
                       Analysis_Type + DTXSID + Species +
                         model + References.Analyzed +
                         Routes + Media +
                         time_units_fitted +
                         AIC ~ param_name,
                       value.var = "Fitted std dev")
  #append "_sd" for these params
  setnames(PK_2comp_sd,
           PK_fit[model %in% "2compartment" &
                    !grepl(x = param_name,
                           pattern = "sigma"), unique(param_name)],
           PK_fit[model %in% "2compartment" &
                    !grepl(x = param_name,
                           pattern = "sigma"), paste0(unique(param_name),
                                                      "_sd")])
  #merge
  PK_2comp <- merge(PK_2comp,
                    PK_2comp_sd,
                    by = intersect(names(PK_2comp),
                                   names(PK_2comp_sd))
  )

  #convert time constants back to hours if necessary
  #Convert units of time constants if rescale_time == TRUE
  time_units <- unique(PK_2comp$time_units_fitted)
  if(!(time_units %in% "hours")){
    time_const <- intersect(names(PK_2comp),
                            c("kelim",
                              "kgutabs",
                              "k12",
                              "k21",
                              "kelim_sd",
                              "kgutabs_sd",
                              "k21_sd",
                              "k12_sd")
    )
    if(time_units == "days"){
      convfun <- function(x) x*24
    }else if(time_units == "weeks"){
      convfun <- function(x) x*(24*7)
    }else if(time_units == "months"){
      convfun <- function(x) x*(24*30)
    }else if(time_units == "years"){
      convfun <- function(x) x*(24*365)
    }else if(time_units == "minutes"){
      convfun <- function(x) x/60
    }
    PK_2comp[, (time_const) := lapply(.SD,
                                      convfun),
             .SDcols = time_const]
  }

  PK_2comp[, time_units_reported := "hours"]

  #in case Fgutabs and V1 were fitted separately, compute Fgutabs/V1
  PK_2comp[is.na(Fgutabs_V1), Fgutabs_V1 := Fgutabs/V1]

  #Total clearance
  PK_2comp[, CLtot := kelim * V1]

  #Total clearance divided by Fgutabs, in case only Fgutabs/V1 was available
  PK_2comp[, CLtot_Fgutabs := kelim / Fgutabs_V1]



  #compute A, B, alpha, beta
  #see https://www.boomer.org/c/p4/c19/c1902.php
  PK_2comp[, alphabeta_sum := kelim + k12 + k21]
  PK_2comp[, alphabeta_prod := kelim * k21]
  PK_2comp[, alpha := (alphabeta_sum + sqrt(alphabeta_sum^2 - 4*alphabeta_prod))/2 ]
  PK_2comp[, beta := (alphabeta_sum - sqrt(alphabeta_sum^2 - 4*alphabeta_prod))/2 ]
  #A and B when only IV data were available (for dose of 1 mg/kg)
  PK_2comp[, A_1mgkg := 1*(alpha - k21) / (V1 * (alpha - beta))]
  PK_2comp[, B_1mgkg:= 1*(k21 - beta) / (V1 * (alpha - beta))]
  #A and B when oral data were available (for dose of 1 mg/kg)
  PK_2comp[!is.na(kgutabs), A_1mgkg := (kgutabs * Fgutabs_V1 *
                     (alpha - k21)) /
             ( (kgutabs - alpha) * (alpha - beta))]
  PK_2comp[!is.na(kgutabs), B_1mgkg := (kgutabs * Fgutabs_V1 *
                     (k21 - beta)) /
             ( (kgutabs - beta) * (alpha - beta))]

  #Apparent volumes of distribution

  #Vbeta
  #Terminal volume of distribution
  #see https://www.boomer.org/c/p4/c19/c1905.php

  PK_2comp[, Vbeta := V1 * kelim / beta]
  PK_2comp[, Vbeta_Fgutabs := (1/Fgutabs_V1) * kelim / beta]

  #Vss
  #apparent volume of distribution at steady state
  #see https://www.boomer.org/c/p4/c19/c1905.php
  PK_2comp[, Vss := V1 * (k21 + k12) / k21]
  #Vss/Fgutabs, in case only Fgutabs/V1 was available and not V1 by itself
  PK_2comp[, Vss_Fgutabs := (1/Fgutabs_V1) * (k21 + k12) / k21]

  #overall relationship: Vbeta > Vss > V1

  #Get Css = average plasma concentration for 1 mg/kg/day every 1 days
  PK_2comp[, Css_iv_1mgkg := 1/(24 * CLtot)]
  PK_2comp[, Css_oral_1mgkg := Fgutabs_V1 * (1/(24*kelim))]

  #half-lives for each phase
  PK_2comp[, halflife_beta := log(2) / beta]
  PK_2comp[, halflife_alpha := log(2) / alpha]
  PK_2comp[, halflife_abs := log(2) / kgutabs]

  #tpeak for oral dose
  #I don't think it can be done analytically
  #so do it numerically
  #use uniroot() to search for zeros of time deriv of Cp vs. time
  #for each set of parameters
  #look between 0 and what would be the 1-compartment tpeak
  #extending the interval if necessary
  #if uniroot() fails, then just return NA
  PK_2comp[!is.na(kgutabs), tpeak_oral := tryCatch(
    uniroot( f = function(x){
    cp_2comp_dt(params = list("kelim" = kelim,
                              "Fgutabs_V1" = Fgutabs_V1,
                              "kgutabs" = kgutabs,
                              "k12" = k12,
                              "k21" = k21),
                time = x,
                dose = 1,
                iv.dose = FALSE)
  },
            lower = 0,
             upper = log(kgutabs / kelim) / (kgutabs - kelim),
  extendInt = "downX", #function should be decreasing
            maxiter = 1000,
  tol = .Machine$double.eps)$root,
  error = function(err) return(NA_real_)),
  by = .(Analysis_Type,
         DTXSID,
         Species,
         model,
         References.Analyzed,
         AIC)]

  #Cpeak for oral dose of 1 mg/kg
  PK_2comp[!is.na(kgutabs),
           Cpeak_oral_1mgkg := cp_2comp(params = list("kelim" = kelim,
                                                  "Fgutabs_V1" = Fgutabs_V1,
                                                  "kgutabs" = kgutabs,
                                                  "k12" = k12,
                                                  "k21" = k21),
                                    time = tpeak_oral,
                                    dose = 1,
                                    iv.dose = FALSE),
           by = .(Analysis_Type,
                  DTXSID,
                  Species,
                  model,
                  References.Analyzed,
                  AIC)]

  #remove temporary calculation columns
  PK_2comp[, c("alphabeta_sum",
               "alphabeta_prod") := NULL]
  PK_out <- copy(PK_2comp)
}

return(PK_out)
}
