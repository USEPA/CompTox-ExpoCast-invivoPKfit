#' Evaluate goodness of fit based on residuals
#'
#' Evaluate multiple metrics of goodness of fit based on residuals
#'
#' For each dataset (chemical and species) and each fitted model, this function
#' first calculates the predictions and residuals, then evaluates the following:
#'
#' - RMSE (root mean squared error)
#' - R-squared of observations vs. predictions
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
#' @param cvt_pre Concentration vs. time data for which to evaluate residuals.
#'   Should be a subset of the table produced by [preprocess_data()], containing
#'   data for the DTXSID and species corresponding to `pk_fit_row`.
#' @param pk_fit_row One row of post-processed fitting output, as produced by
#'   [merge_fits()]. This corresponds to one DTXSID, species, analysis type, and
#'   one set of fitting options.
#' @return A `data.table` of goodness-of-fit measures for each dataset, each
#'   model, and each analysis (joint, separate, and pooled)
#' @author Caroline Ring
#' @export
#' @import data.table
evaluate_resids <- function(cvt_pre,
                            pk_fit_row,
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
cvt_pre[, pred_conc := predict_conc(cvt_pre = .SD,
                               pk_fit_row = pk_fit_row)]

 cvt_pre[, pred_conc_Dose := pred_conc/Dose]


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

#Compute residuals in different ways:
 #log-transformed yes/no, dose-normalized yes/no
 #Do this regardless of the original fitting conditions,
 #so we can compare residuals across fitting conditions.
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
    #for each group, do this:
    #Create a list of arguments to calc_rmse() and calc_rsq()
    #This will consist of different sets of columns in pred_DT,
    #depending on whether dose_norm is TRUE or FALSE
    if(dose_norm %in% FALSE){
      #Grab non-dose-normalized columns to pass to calc_rmse() and calc_rsq()
      tmp <- list("group_mean" = Conc,
                  "group_sd" = Conc_SD,
                  "group_n" = N_Subjects,
                  "group_LOQ" = LOQ,
                  "pred" = pred_conc)

      #Calculate prediction normalized by concentration
      #this will give you "within a factor of two" info
      if(log_trans %in% FALSE){
      pred_obs <- pred_conc/Conc
      }else{
        pred_obs <- log(pred_conc)/log(Conc)
      }
    }else{
      #Grab dose-normalized columns to pass to calc_rmse() and calc_rsq()
      tmp <- list("group_mean" = Conc_Dose,
                  "group_sd" = Conc_SD_Dose,
                  "group_n" = N_Subjects,
                  "group_LOQ" = LOQ_Dose,
                  "pred" = pred_conc_Dose)
      #Calculate prediction normalized by concentration
      #this will give you "within a factor of two" info
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
#create an output table: RMSE, Rsq, and the mean, median, min, and max of
#pred/obs
    outDT <- data.table("RMSE" = RMSE,
         "rsq" = rsq,
         "pred_obs_mean" = mean(pred_obs),
         "pred_obs_median" = median(pred_obs),
         "pred_obs_min" = min(pred_obs),
         "pred_obs_max" =max(pred_obs))

  #To the resulting data.table, add variables recording dose_norm and log_trans
  outDT[, dose_norm := dose_norm]
  outDT[, log_trans := log_trans]
  outDT
},
SIMPLIFY = FALSE)

#the result is a list of data.tables, one for each type of residual
#rowbind them together
resid_DT <- rbindlist(resid_DT_list,
                      use.names = TRUE,
                      fill = TRUE)

 #How big are "overall" error SDs?

 #Get overall error SDs (pseudo-RMSE?) for various analysis types & models.
 #This is the best-fit standard deviation for the log-likelihood, evaluated
 #holding the PK model parameters constant.

 #Again, use mapply() to loop over combinations of log_trans and dose_norm
 error_sd_DT_list <- mapply(function(log_trans,
                    dose_norm){
   if(verbose %in% TRUE){
     message(paste0("Estimating overall error sigma for:\n",
                    "log transformation ", log_trans, "\n",
                    "dose normalization ", dose_norm, "\n"))
   }

    #get "pooled" error SD estimate keeping params fixed at their fitted values
    #do this for all four combos of log TRUE/FALSE and dose-normalization
    #TRUE/FALSE
    out_DT <- get_error_sd(pk_fit_row = pk_fit_row,
                 newdata = newdata,
                 model = model,
                 dose_norm = dose_norm,
                 fit_log_conc = log_trans,
                 optimx_args = optimx_args)

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
                   by = c("dose_norm",
                          "log_trans"))

return(eval_DT)

}

