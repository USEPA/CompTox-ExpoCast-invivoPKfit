#' Evaluate goodness of fit
#'
#' Compute multiple metrics to evaluate how well the various TK models fit the
#' data.
#'
#' @param cvt_pre Original concentration vs. time data, pre-processed using
#'   `preprocess_data()`.
#' @param fit_flat Fitting output for the flat model: the output of `fit_all()`
#'   with `model = "flat"`
#' @param fit_1comp Fitting output for the 1-compartment model: the output of
#'   `fit_all()` with `model = "1compartment"`
#' @param fit_2comp Fitting output for the 2-compartment model: the output of
#'   `fit_all()` with `model = "2compartment"`
#' @param fit_dose_norm Logical: Whether the model fits were done on
#'   dose-normalized data (TRUE), or not (FALSE). This needs to be the same as
#'   whatever you used when calling `fit_all()`, or the goodness-of-fit metrics
#'   will be incorrect!
#' @return A `data.table` of goodness-of-fit measures for each dataset, each
#'   model, and each analysis (joint, separate, and pooled)
#' @author Caroline Ring
#' @export
#' @import data.table
evaluate_fit <- function(cvt_pre,
                         fit_flat,
                         fit_1comp,
                         fit_2comp,
                         fit_dose_norm){




  post_flat <- postprocess_data(PK_fit = fit_flat,
                                model = "flat")
  post_1comp <- postprocess_data(PK_fit = fit_1comp,
                                 model = "1compartment")
  post_2comp <- postprocess_data(PK_fit = fit_2comp,
                                 model = "2compartment")

#Merge together the post-processed data, using the following set of key columns
  idcols <- c( "Analysis_Type",
               "DTXSID",
               "Species",
               "References.Analyzed",
               "Studies.Analyzed",
               "N_Routes",
               "N_Media",
               "time_units_reported")
#Change the names of the parameter columns to include the model: e.g.
#"kelim.1compartment", "kelim.2compartment", etc.
  setnames(post_1comp,
           setdiff(names(post_1comp),
                   idcols),
           paste(setdiff(names(post_1comp),
                         idcols),
                 "1compartment",
                 sep = "."))

  setnames(post_2comp,
           setdiff(names(post_2comp),
                   idcols),
           paste(setdiff(names(post_2comp),
                         idcols),
                 "2compartment",
                 sep = "."))

  setnames(post_flat,
           setdiff(names(post_flat),
                   idcols),
           paste(setdiff(names(post_flat),
                         idcols),
                 "flat",
                 sep = "."))
#Merge together all the tables of fitted parameters
  pk_fit <- merge(post_2comp,
                  merge(post_flat, post_1comp, by = idcols),
                  by = idcols)

  #find the winning model for each dataset (minimum AIC)
  pk_AIC <- melt(pk_fit,
                 id.vars = idcols,
                 measure.vars = c("AIC.flat",
                                  "AIC.1compartment",
                                  "AIC.2compartment"),
                 value.name = "AIC",
                 variable.name = "model")

  pk_AIC[, model := gsub(x = model,
                         pattern = "AIC.",
                         replacement = "",
                         fixed = TRUE)]


  #sort AIC from smallest to largest within each dataset, putting NAs last
  setorderv(pk_AIC,
            c(idcols, "AIC"),
            na.last = TRUE)
  #the winning model will now be listed first within each dataset (smallest AIC)
  model_win <- pk_AIC[,
                      .(winmodel = model[1]),
                      by = idcols]

  #Merge in the table of winning models
  pk_fit <- model_win[pk_fit,
                      on = idcols]

  #Merge the table of fitted parameters with the original CVT data.

  #Ensure comma-separated lists of studies analyzed are properly sorted
  pk_fit[, Studies.Analyzed := {
    tmp <- strsplit(x = Studies.Analyzed,
                    split = ", ")
    sapply(tmp,
           function(x) paste(sort(x), collapse = ", "))
  }]

  #Create a column with pooled studies in the original CVT data: comma-separated
  #list of all studies for a given DTXSID and Species
  cvt_pre[, Study_pooled :=
        paste(
          sort(unique(Study)),
          collapse = ", "
        ),
      by = .(DTXSID, Species)]

  #Melt to longer format. The effect will be as though we row-bound two versions
  #of DF together: one with the original single references, and one with
  #comma-separated combined References. There is one column "Study_type" that
  #contains either "Study" or "Study_pooled", and then one column "Study" that
  #contains the studyID or the comma-separated list of study IDs.
  cvt_pre[, Study_orig := Study]
  cvt1 <- melt(cvt_pre,
               measure.vars = c("Study", "Study_pooled"),
               variable.name = "Study_type",
               value.name = "Studies.Analyzed")

  #There will be duplicated rows for single-study analyses, because
  #"Study_pooled" (the comma-separated list of study IDs) is the same as "Study"
  #(the individual study ID) for these analyses. Remove the duplicated rows.
  cvt1 <- unique(cvt1[, .(DTXSID, Compound, Species, Studies.Analyzed, Reference,
                          Route, Dose, iv, Media, Study_orig, Time,
                          Value, LOQ, Value_SD, N_Subjects,
                          Conc, Conc_Dose, Detect, Value_Dose,
                          Value_SD_Dose, LOQ_Dose)])

  #merge original CVT data with fitted PK parameters
  cvt_fit <- pk_fit[cvt1,
                       on = c("DTXSID" = "DTXSID",
                              "Species" = "Species",
                              "Studies.Analyzed" = "Studies.Analyzed"),
                       allow.cartesian = TRUE]

  #Also, prepare a table of study-specific sigmas

  #These are from the original full tables of fitting outputs
  idcols <- c( "Analysis_Type",
               "DTXSID",
               "Species",
               "References.Analyzed",
               "Studies.Analyzed")

  fit_all <- rbindlist(list(fit_1comp, fit_2comp, fit_flat),
                      use.names = TRUE,
                      fill = TRUE)
  fit_all <- fit_all[grepl(x = param_name,
                         pattern = "sigma"),
                   .SD, .SDcols= c(idcols, "model",
                                   "param_name",
                                   "Fitted mean",
                                   "Fitted std dev")]
  fit_melt <- dcast(fit_all,
                   formula = Analysis_Type + DTXSID + Species +
                     References.Analyzed +
                     Studies.Analyzed + param_name ~ model,
                   value.var = "Fitted mean")
  setnames(fit_melt,
           unique(fit_all$model),
           paste0("sigma.", unique(fit_all$model)))
  fit_melt[param_name %in% "sigma", sigma_studyID := Studies.Analyzed]
  fit_melt[!(param_name %in% "sigma"),
           sigma_studyID := gsub(x = param_name,
                                 pattern = "sigma_study_",
                                 replacement = "",
                                 fixed = TRUE)]
  fit_melt_sd <- dcast(fit_all,
                      formula = Analysis_Type + DTXSID + Species +
                        References.Analyzed + Studies.Analyzed +
                        param_name ~ model,
                      value.var = "Fitted std dev")
  setnames(fit_melt_sd,
           unique(fit_all$model),
           paste0("sigma_sd.", unique(fit_all$model)))
  fit_melt_sd[param_name %in% "sigma", sigma_studyID := Studies.Analyzed]
  fit_melt_sd[!(param_name %in% "sigma"),
              sigma_studyID := gsub(x = param_name,
                                    pattern = "sigma_study_",
                                    replacement = "",
                                    fixed = TRUE)]

  fit_melt_all <- merge(fit_melt,
                        fit_melt_sd,
                        by = intersect(names(fit_melt),
                                                           names(fit_melt_sd)))

  #Get SD of data for each study
 cvt_fit[, sd_Value_Dose := sd(Value_Dose, na.rm = TRUE),
                    by = Study_orig]

  #Merge SD of data with fitted SDs
  cvt_fit <- cvt_fit[fit_melt_all,
                     on = c("Study_orig" = "sigma_studyID")]

  #Now that we have all the data in order, start calculating goodness of fit metrics

  #Model predicted values for each model
  for(this_model in c("flat", "1compartment", "2compartment")){
    predname <- paste("pred",
                      this_model,
                      sep = ".")
    modelfun <- ifelse(this_model %in% "flat",
                       cp_flat,
                       ifelse(this_model %in% "1compartment",
                              cp_1comp,
                              cp_2comp))
    #get the model parameter names for this model
    param_cols <- get_model_paramnames(model = this_model)
    param_cols <- paste(param_cols,
                        this_model,
                        sep = ".")
    #Evaluate model predictions for this model for each dataset & analysis
    cvt_fit[, (predname) :=  do.call(what = modelfun,
                                     args = c(.SD,
                                              list(time = Time,
                                                   dose = Dose,
                                                   iv.dose = Route %in% "iv",
                                                   medium = Media)) ),
            .SDcols = param_cols,
            by = .(DTXSID,
                   Species,
                   Analysis_Type)]
  }

  #Melt to longer format by model
  cvt_melt <- melt(cvt_fit,
                   id.vars = idcols,
                   measure.vars = c("pred.flat",
                                    "pred.1compartment",
                                    "pred.2compartment"),
                   variable.name = "model",
                   value.name = "pred")
  cvt_melt[, model := gsub(x = model,
                           pattern = "pred.",
                           replacement = "",
                           fixed = TRUE)]

    #calculate dose-normalized predictions
    cvt_fit[, "pred_Dose" :=  pred / Dose ]

    #Get residuals
   #Non-dose-normalized
    cvt_fit[, resid := Value - pred]
    #Dose-normalized
    cvt_fit[, resid_Dose :=  Value_Dose - pred_Dose ]

    # Root mean squared error
    ## Non-dose-normalized
    ### Calculate RMSE per chemical and species (i.e., per dataset)
    cvt_fit[, RMSE_chemical_species := sqrt(mean(resid^2)),
            by = .(DTXSID, Species, Analysis_Type, model)]

    #### Calculate RMSE per dataset & study
    cvt_fit[, RMSE_chemical_species_study := sqrt(mean(resid^2)),
            by = .(DTXSID, Species, Study, Analysis_Type, model)]

    ### Calculate RMSE per dataset & dose
    cvt_fit[, RMSE_chemical_species_dose := sqrt(mean(resid^2)),
            by = .(DTXSID, Species, Dose, Analysis_Type, model)]

 #### Calculate RMSE per dataset, study, and dose
    cvt_fit[, RMSE_chemical_species_study_dose := sqrt(mean(resid^2)),
            by = .(DTXSID, Species, Study, Dose, Analysis_Type, model)]

    ## Dose-normalized RMSE
    ### Calculate RMSE per chemical and species (i.e., per dataset)
    cvt_fit[, RMSE_Dose_chemical_species := sqrt(mean(resid_Dose^2)),
            by = .(DTXSID, Species, Analysis_Type, model)]

    #### Calculate RMSE per dataset & study
    cvt_fit[, RMSE_Dose_chemical_species_study := sqrt(mean(resid_Dose^2)),
            by = .(DTXSID, Species, Study, Analysis_Type, model)]

    #### Calculate RMSE per dataset & dose
    cvt_fit[, RMSE_Dose_chemical_species_dose := sqrt(mean(resid_Dose^2)),
            by = .(DTXSID, Species, Dose, Analysis_Type, model)]

    #### Calculate RMSE per dataset, study, and dose
    cvt_fit[, RMSE_Dose_chemical_species_study_dose := sqrt(mean(resid_Dose^2)),
            by = .(DTXSID, Species, Study, Dose, Analysis_Type, model)]


    # R-squared of lm(Observed ~ Predicted)

    ## Non-dose-normalized
    ### Calculate R^2 of fit per dataset
    cvt_fit[, R2_chemical_species := summary(lm(Value ~ pred))$r.squared,
            by = .(DTXSID, Species, Analysis_Type, model)]

    ### Calculate R^2 of fit per dataset & study
    cvt_fit[, R2_chemical_species_study := summary(lm(Value ~ pred))$r.squared,
            by = .(DTXSID, Species, Study, Analysis_Type, model)]

    ### Calculate R^2 of fit per dataset & dose
    cvt_fit[, R2_chemical_species_dose := summary(lm(Value ~ pred))$r.squared,
            by = .(DTXSID, Species, Dose, Analysis_Type, model)]

    ### Calculate R^2 of fit per dataset, study,  & dose
    cvt_fit[, R2_chemical_species_study_dose := summary(lm(Value ~ pred))$r.squared,
            by = .(DTXSID, Species, Study, Dose, Analysis_Type, model)]

    ## Dose-normalized
    ### Calculate R^2 of fit per dataset
    cvt_fit[, R2_Dose_chemical_species := summary(lm(Value_Dose ~ pred_Dose))$r.squared,
            by = .(DTXSID, Species, Analysis_Type, model)]

    ### Calculate R^2 of fit per dataset & study
    cvt_fit[, R2_Dose_chemical_species_study := summary(lm(Value_Dose ~ pred_Dose))$r.squared,
            by = .(DTXSID, Species, Study, Analysis_Type, model)]

    ### Calculate R^2 of fit per dataset & dose
    cvt_fit[, R2_Dose_chemical_species_dose := summary(lm(Value_Dose ~ pred_Dose))$r.squared,
            by = .(DTXSID, Species, Dose, Analysis_Type, model)]

    ### Calculate R^2 of fit per dataset, study,  & dose
    cvt_fit[, R2_Dose_chemical_species_study_dose := summary(lm(Value_Dose ~ pred_Dose))$r.squared,
            by = .(DTXSID, Species, Study, Dose, Analysis_Type, model)]

  #Calculate the empirical Cpeak/tpeak
  #Non-dose-normalized
  cvt_fit[!is.na(Value),
          c("tpeak", "Cpeak") := get_peak(x = Time,
                                            y = Value),
          by = .(DTXSID, Species, Analysis_Type, model)]
  #Dose-normalized
  cvt_fit[!is.na(Value),
          c("tpeak_Dose", "Cpeak_Dose") := get_peak(x = Time,
                                            y = Value_Dose),
          by = .(DTXSID, Species, Analysis_Type, model)]

  #Calculate

  #



}
