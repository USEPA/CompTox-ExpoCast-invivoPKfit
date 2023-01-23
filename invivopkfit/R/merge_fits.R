merge_fits <- function(fit_flat,
                       fit_1comp,
                       fit_2comp){

  fitlist <- list("flat" = fit_flat,
                  "1compartment" = fit_1comp,
                  "2compartment" = fit_2comp)
  pk_fit <- rbindlist(lapply(seq_along(fitlist),
                     function(i) {
                       this_model <- names(fitlist)[i]
                       this_post <- postprocess_data(PK_fit = fitlist[[i]],
                                                  model = this_model)
                       return(this_post)
                     }
  ),
  use.names = TRUE,
  fill = TRUE)

  #Merge together the post-processed data, using the following set of key columns
  idcols <- c( "Analysis_Type",
               "DTXSID",
               "Species",
               "References.Analyzed",
               "Studies.Analyzed",
               "N_Routes",
               "N_Media",
               "time_units_reported",
               "rescale_time",
               "fit_conc_dose",
               "fit_log_conc")

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
