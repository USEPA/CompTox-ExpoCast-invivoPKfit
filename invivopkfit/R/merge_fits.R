merge_fits <- function(fit_flat,
                       fit_1comp,
                       fit_2comp){

  #row-bind together the fit tables from each model
  fitlist <- list("flat" = fit_flat,
                  "1compartment" = fit_1comp,
                  "2compartment" = fit_2comp)

  #post-process each fit table
  postlist <-  lapply( #post-process each fit table before row-binding
    seq_along(fitlist),
    function(i) {
      this_model <- names(fitlist)[i]
      this_post <- postprocess_data(
        PK_fit = fitlist[[i]],
        model = this_model
      )
      return(this_post)
    }
  )

  names(postlist) <- names(fitlist)

  #now row-bind post-processed tables
  pk_fit <- data.table::rbindlist(postlist,
    use.names = TRUE,
    fill = TRUE)

#Define a set of "ID columns" that define a specific fit (for a single dataset,
#analysis type, model, & fitting conditions)
  idcols <- c("Analysis_Type",
               "DTXSID",
               "Species",
               "References.Analyzed",
               "Studies.Analyzed",
               "N_Routes",
               "N_Media",
               "time_units_reported",
               "rescale_time",
               "fit_conc_dose",
               "fit_log_conc",
               "method")

  #Calculate "p-value" for 1-comp and 2-comp models vs. null model
  #First get table of AICs by dataset, analysis, model
  fit_all <- rbindlist(fit_list,
                       use.names = TRUE,
                       fill = TRUE)
  AIC_DT <- unique(fit_all[, .SD, .SDcols = c(idcols,
                                              "model",
                                              "AIC")])
  #cast wider to get AIC by model for each analysis
  AIC_DT_wide <- dcast(AIC_DT,
                       ... ~ model,
                       value.var = "AIC")
  #calculate p-value
  AIC_DT_wide[, pval_1compartment := exp((`1compartment` - flat)/2)]
  AIC_DT_wide[, pval_2compartment := exp((`2compartment` - flat)/2)]
  AIC_DT_wide[, pval_flat := 1]
  #melt p-values longer again
  pval_DT <- melt(AIC_DT_wide,
                  id.vars = idcols,
                  measure.vars = c("pval_1compartment",
                                   "pval_2compartment",
                                   "pval_flat"),
                  variable.name = "model",
                  value.name = "AIC_pval_vs_flat",
                  variable.factor = FALSE)
  #identify by model
  pval_DT[, model := gsub(x =model,
                          pattern = "pval_",
                          replacement = "",
                          fixed = TRUE)]
  #merge in
  pk_fit <- merge(pk_fit,
                  pval_DT,
                  by = c(idcols, "model"))

  #find the winning model for each dataset (minimum AIC)

  #sort AIC from smallest to largest within each dataset, putting NAs last
  setorderv(pk_fit,
            c(idcols, "AIC"),
            na.last = TRUE)

  #the winning model will now be listed first within each dataset (smallest AIC)
  #add a column noting the winning model
  pk_fit[, winmodel := model[1], by = idcols]
  #In case no model fit was done (e.g. no detected data), note that fact
  pk_fit[, no_fit := all(is.na(AIC)), by = idcols]
  pk_fit[no_fit %in% TRUE, winmodel := "None (no fit)"]
  pk_fit[, no_fit := NULL]

 #Cast to wide format
  #Produce a formula to use for casting
  #[idcols] ~ model
  #add "winmodel" to idcols
  idcols <- c(idcols, "winmodel")
  cast_formula <- as.formula(
    paste(
      paste(idcols, collapse = " + "),
      "model",
      sep = " ~ "
    )
  )

  #column names to be used as value variables
  #these will come out named as [value_name].[model]
  value_names <- setdiff(names(pk_fit),
                         c(idcols,
                           "model"))
#do the cast
  pk_fit <- dcast(pk_fit,
                  formula = as.formula(cast_formula),
                  value.var = value_names,
                  sep = ".",
                  drop = TRUE)

  #drop any columns that are all-NA
  #these are parameter & model combinations that do not occur, like "V1.flat"
  na_cols <- pk_fit[, sapply(.SD, function(x) all(is.na(x)))]
  na_colnames <- names(pk_fit)[na_cols]
  pk_fit[, (na_colnames) := NULL]

  return(pk_fit)
}
