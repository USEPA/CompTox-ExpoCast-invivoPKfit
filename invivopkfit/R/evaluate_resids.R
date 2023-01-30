#' Evaluate goodness of fit based on residuals
#'
#' Evaluate multiple metrics of goodness of fit based on residuals
#'
#' For each dataset (chemical and species) and each fitted model, this function
#' first calculates the predictions and residuals, then evaluates the following:
#'
#' - RMSE (root mean squared error)
#' - R-squared of observations vs. predictions
#' - The p-value from the Breusch-Pagan test of heteroscedasticity
#' - Estimated "overall" error standard deviation (details below)
#'
#' Unique fitted models are defined by unique combination of analysis type
#' ("Joint", "Separate", or "Pooled"), studies analyzed, and model ("flat",
#' "1compartment", or "2compartment").
#'
#' # Root mean squared error (RMSE)
#'
#' RMSE is calculated in [calc_rmse()].
#'
#' # Estimated overall error standard deviation
#'
#' The "overall" error standard deviation is the value for residual error
#' standard deviation ("sigma") that minimizes the log-likelihood function when
#' the TK model parameters are fixed at their fitted values. This is essentially
#' a more rigorous RMSE, because it correctly models the likelihood of
#' non-detect (below-LOQ) observations as the cumulative probability that the
#' model prediction is below the reported LOQ, rather than making any
#' substitutions for non-detects.
#'
#' @param cvt_pre Original concentration vs. time data, pre-processed using
#'   `preprocess_data()`.
#' @param fit_flat Fitting output for the flat model: the output of `fit_all()`
#'   with `model = "flat"`
#' @param fit_1comp Fitting output for the 1-compartment model: the output of
#'   `fit_all()` with `model = "1compartment"`
#' @param fit_2comp Fitting output for the 2-compartment model: the output of
#'   `fit_all()` with `model = "2compartment"`
#' @param fit_conc_dose Logical: Whether the model fits were done on
#'   dose-normalized data (TRUE), or not (FALSE). This needs to be the same as
#'   whatever you used when calling `fit_all()`, or the goodness-of-fit metrics
#'   will be incorrect!
#' @param fit_log_conc Logical: Whether the model fits were done on
#'   log-transformed data (TRUE), or not (FALSE). This needs to be the same as
#'   whatever you used when calling `fit_all()`, or the goodness-of-fit metrics
#'   will be incorrect!
#' @return A `data.table` of goodness-of-fit measures for each dataset, each
#'   model, and each analysis (joint, separate, and pooled)
#' @author Caroline Ring
#' @export
#' @import data.table
evaluate_resids <- function(cvt_pre,
                         fit_flat,
                         fit_1comp,
                         fit_2comp,
                         fit_conc_dose = TRUE,
                         fit_log_conc = FALSE,
                         optimx_args = list(
                           "method" = "bobyqa",
                           "itnmax" = 1e6,
                           "control" = list("kkt" = FALSE)
                         ),
                         verbose = TRUE){


#Calculate RMSE, R-squared of obs vs pred,
  #and heteroscedasticity of residuals

#First: calculate residuals

#Get predicted concentrations & AUCs
  if(verbose %in% TRUE){
    message("Getting predicted concentrations...")
  }
 pred_DT <- get_all_predictions(cvt_pre =cvt_pre,
                               fit_flat = fit_flat,
                               fit_1comp = fit_1comp,
                               fit_2comp = fit_2comp)

 pred_DT[, pred_conc_Dose := pred_conc/Dose]

 pred_DT[, n_obs := .N, by = .(DTXSID, Species,
                               Analysis_Type,
                               Studies.Analyzed,
                               model)]

 pred_DT[, n_detect := sum(Detect %in% "Detect"),
         by = .(DTXSID, Species,
                               Analysis_Type,
                               Studies.Analyzed,
                               model)]

#For the diagnostics: We need to do all of these under all four conditions on
#log-scaling TRUE/FALSE and dose-normalization TRUE/FALSE, no matter how the
#model itself was fit. That way, we can compare apples to apples across fitting
#methods.

#I think that the non-log, non-dose-normalized metrics can be a "baseline" but
#I'm not sure. Better to calculate them all.

#This is somewhat time-consuming but I think it will be OK.
#Calculate RMSE, Rsquared
 #Skip B-P test of heteroskedasticity for now
 if(verbose %in% TRUE){
   message("Calculating RMSE, r-squared...")
 }

 #Make a table of conditions on log transformation and dose normalization
resid_cond <- as.data.table(expand.grid(log_trans = c(TRUE, FALSE),
                          dose_norm = c(TRUE, FALSE)))

#use mapply() to loop over these combinations of conditions
resid_DT_list <- mapply(function(log_trans,
                dose_norm){
  if(verbose %in% TRUE){
    message(paste0("Calculating RMSE and r-squared for:\n",
                   "log transformation ", log_trans, "\n",
                   "dose normalization ", dose_norm, "\n"))
  }
  #Use data.table syntax to loop over groups within pred_DT:
  #unique combinations of DTXSID, Species,
  #Analysis_Type, Studies.Analyzed,
  #model (see "by =" argument)
  outDT <- pred_DT[, {
    #for each group, do this:
    #Create a list of arguments to calc_rmse() and calc_rsq()
    #This will consist of different sets of columns in pred_DT,
    #depending on whether dose_norm is TRUE or FALSE
    if(dose_norm %in% FALSE){
      #Grab non-dose-normalized columns
      tmp <- list("group_mean" = Conc,
                  "group_sd" = Conc_SD,
                  "group_n" = N_Subjects,
                  "group_LOQ" = LOQ,
                  "pred" = pred_conc)
      if(log_trans %in% FALSE){
      pred_obs <- pred_conc/Conc
      }else{
        pred_obs <- log(pred_conc)/log(Conc)
      }
    }else{
      #Grab dose-normalized columns
      tmp <- list("group_mean" = Conc_Dose,
                  "group_sd" = Conc_SD_Dose,
                  "group_n" = N_Subjects,
                  "group_LOQ" = LOQ_Dose,
                  "pred" = pred_conc_Dose)
      if(log_trans %in% FALSE){
       pred_obs <- pred_conc_Dose/Conc_Dose
      }else{
        pred_obs <- log(pred_conc_Dose)/log(Conc_Dose)
      }
    }
    #if both observed and predicted < LOQ, set pred_obs to 1
    pred_obs[Conc <= LOQ &
                 pred_conc <= LOQ] <- 1

    #Calculate RMSE and rsq
      RMSE <- do.call(calc_rmse,
                      args = c(tmp,
                               "log" = log_trans))
      rsq <- do.call(calc_rsq,
                     args = c(tmp,
                              "log" = log_trans))

    list("RMSE" = RMSE,
         "rsq" = rsq,
         "pred_obs_mean" = mean(pred_obs),
         "pred_obs_min" = min(pred_obs),
         "pred_obs_max" =max(pred_obs))

  }, by = .(DTXSID, Species,
         Analysis_Type, Studies.Analyzed, model,
         n_obs,
         n_detect)]

  #To the resulting data.table, add variables recording dose_norm and log_trans
  outDT[, dose_norm := dose_norm]
  outDT[, log_trans := log_trans]
  return(outDT)
},
dose_norm = resid_cond$dose_norm,
log_trans = resid_cond$log_trans,
SIMPLIFY = FALSE)

#the result is a list of data.tables
#rowbind them together
resid_DT <- rbindlist(resid_DT_list,
                      use.names = TRUE,
                      fill = TRUE)

 #How big are "overall" error SDs?

 #Get overall error SDs (pseudo-RMSE?) for various analysis types & models.
 #This is the best-fit standard deviation for the log-likelihood, evaluated
 #holding the PK model parameters constant.

#merge together all tables of fitted parameters
pk_fit <- merge_fits(fit_flat = fit_flat,
                     fit_1comp = fit_1comp,
                     fit_2comp = fit_2comp)

 #Enumerate combinations of chemical, species, analysis type, studies analyzed
 fit_combs <- unique(fit_flat[, .(DTXSID,
                                 Species,
                                 Analysis_Type,
                                 Studies.Analyzed)])
#add models to the table
 fit_combs <- fit_combs[, .(model = c("flat",
                                      "1compartment",
                                      "2compartment")),
                        by = names(fit_combs)]

 setnames(fit_combs,
          names(fit_combs),
          paste(names(fit_combs),
                "in",
                sep = "_"))

 #Again, use mapply() to loop over combinations of log_trans and dose_norm
 error_sd_DT_list <- mapply(function(log_trans,
                    dose_norm){
   if(verbose %in% TRUE){
     message(paste0("Estimating overall error sigma for:\n",
                    "log transformation ", log_trans, "\n",
                    "dose normalization ", dose_norm, "\n"))
   }
     out_DT <- fit_combs[, {
    studies_list <- strsplit(Studies.Analyzed_in,
                             split = ", ")[[1]]
    #get data.table of newdata to be predicted
    #leave out DTXSID, Species columns
    #since these will be added back in from fit_combs
    newdata <- cvt_pre[DTXSID %in% DTXSID_in &
                         Species %in% Species_in &
                         Study %in% studies_list,
                       .SD,
                       .SDcols = setdiff(names(cvt_pre),
                                         c("DTXSID", "Species"))]
    #get "pooled" error SD estimate keeping params fixed at their fitted values
    #do this for all four combos of log TRUE/FALSE and dose-normalization
    #TRUE/FALSE
    get_error_sd(pk_fit = pk_fit,
                 newdata = newdata,
                 DTXSID_in = DTXSID_in,
                 Species_in = Species_in,
                 Analysis_Type_in = Analysis_Type_in,
                 Studies.Analyzed_in = Studies.Analyzed_in,
                 model_in = model_in,
                 fit_conc_dose = dose_norm,
                 fit_log_conc = log_trans,
                 optimx_args = optimx_args)
  },
  by = .(DTXSID_in,
         Species_in,
         Analysis_Type_in,
         Studies.Analyzed_in,
         model_in)]

  #remove the "_in" from variable names
  setnames(out_DT,
           names(out_DT),
           gsub(x = names(out_DT),
                pattern = "_in",
                replacement = ""))

  #Add columns recording dose_norm and log_trans
  out_DT[, dose_norm := dose_norm]
  out_DT[, log_trans := log_trans]
  return(out_DT)
 },
 dose_norm = resid_cond$dose_norm,
 log_trans = resid_cond$log_trans,
 SIMPLIFY = FALSE)

 #the result is a list of data.tables --
 #rowbind them together
 error_sd_DT <- rbindlist(error_sd_DT_list,
                          use.names = TRUE,
                          fill = TRUE)

#merge together the diagnostics
  eval_DT <- merge(resid_DT, error_sd_DT,
                   by = c("DTXSID",
                          "Species",
                          "Analysis_Type",
                          "Studies.Analyzed",
                          "model",
                          "dose_norm",
                          "log_trans"))

  if(verbose %in% TRUE){
    message("Getting AICs and log-likelihoods...")
  }
  fit_info <- rbindlist(
    list(
      unique(fit_flat[, .(DTXSID,
                                         Species,
                                         Analysis_Type,
                                         Studies.Analyzed,
                                         model,
                                         AIC,
                                         LogLikelihood)]),
                       unique(fit_1comp[, .(DTXSID,
                                           Species,
                                           Analysis_Type,
                                           Studies.Analyzed,
                                           model,
                                           AIC,
                                           LogLikelihood)]),
                       unique(fit_2comp[, .(DTXSID,
                                           Species,
                                           Analysis_Type,
                                           Studies.Analyzed,
                                           model,
                                           AIC,
                                           LogLikelihood)])
      ),
    use.names = TRUE,
    fill = TRUE
  )

  eval_DT <- merge(fit_info,
                   eval_DT,
                   by = c("DTXSID",
                          "Species",
                          "Analysis_Type",
                          "Studies.Analyzed",
                          "model"))


return(eval_DT)

}

