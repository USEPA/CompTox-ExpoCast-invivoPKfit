#' ---
#' title: "Evaluate fits"
#' author: "Caroline Ring"
#' date: "2023-01-23"
#' output: html_document
#' ---
#' 
## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(results = TRUE)
knitr::opts_chunk$set(message = FALSE)

#' 
## ---- include = FALSE----------------------------------------------------
library(tidyverse)
library(data.table)
library(readxl)
devtools::load_all("../invivopkfit")

#' 
#' # Introduction
#' 
#' This document contains code to evaluate the goodness-of-fit of the models fit to data in `MainAnalysis.Rmd`.
#' 
#' # Import and preprocess CvTdb data
#' 
#' ## Import CvTdb data
#' 
#' Read in CvTdb data that was pulled and normalized in `pulling_oral_iv_data_cvtdb.Rmd`.
#' 
## ---- load_cvt, eval = TRUE----------------------------------------------
cvt_data <- fread("../inst/ext/cvt_data_2023_01_05.csv")

#' 
#' ## Timestamp
#' Timestamp the file name with the current date and time (at the beginning of the runs) in format YYYY_MM_DD_hh_mm.
#' 
## ------------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y_%m_%d")

#' 
#' ## Settings for pre-processing data
#' 
#' List of new_name = old_name to supply for preprocessing data
#' 
## ------------------------------------------------------------------------
names_list <- list(
  "Compound_Dosed" = "studies.test_substance_name_original",
  "DTXSID_Dosed" = "chemicals_dosed.dsstox_substance_id",
  "CAS_Dosed" = "chemicals_dosed.dsstox_casrn",
  "Compound_Analyzed" = "series.analyte_name_original",
  "DTXSID_Analyzed" = "chemicals_analyzed.dsstox_substance_id",
  "CAS_Analyzed" = "chemicals_analyzed.dsstox_casrn",
  "Reference" = "documents_reference.id",
  "Extraction" = "documents_extraction.id",
  "Species" = "subjects.species",
  "Weight" ="subjects.weight_kg",
  "Weight.Units" = NULL,
  "Dose" = "studies.dose_level_normalized_corrected", #note that dose levels were corrected in pulling_oral_iv_data_cvtdb.Rmd
  "Dose.Units" = NULL,
  "Time" = "conc_time_values.time_hr",
  "Time.Units" = NULL,
  "Media" = "series.conc_medium_normalized",
  "Value" = "conc_time_values.conc", #already normalized to mg/L units in CvTdb
  "Value.Units" = NULL,
  "Route" = "studies.administration_route_normalized",
  "LOQ" = "series.loq_normalized", #note that LOQs were normalized to mg/L units in pulling_oral_iv_data_cvtdb.Rmd
  "Subject" = "subjects.id",
  "N_Subjects" = "series.n_subjects_in_series",
  "Value_SD" = "conc_time_values.conc_sd_normalized", #note that SDs were normalized to mg/L units in pulling_oral_iv_data_cvtdb.Rmd
  "Study_ID" = "studies.id",
  "Series_ID" = "series.id"
  )

#' 
#' List of default values for new columns added in preprocessing data
#' 
## ------------------------------------------------------------------------
defaults_list <- list(
                      "Weight.Units" = "kg",
                      "Dose.Units" = "mg/kg",
                      "Time.Units" = "hours",
                      "Value.Units" = "mg/L")

#' 
#' ## Preprocess data
#' 
## ------------------------------------------------------------------------
cvt <- preprocess_data(data.set = cvt_data,
                       names_list =names_list,
                       defaults_list =   defaults_list,
                       ratio_conc_to_dose = 1,
                       calc_loq_factor = 0.45,
                       routes_keep = c("po", "iv"),
                       media_keep = c("blood", "plasma"),
                       impute_loq = TRUE,
                       impute_sd = TRUE,
                       study_def = c("DTXSID", "Species", "Reference", "Route", "Media"),
                       suppress.messages = TRUE)

#save pre-processed data
# write.csv(cvt,
#           paste0("../inst/ext/cvt_preprocessed_",
#                  timestamp,
#           ".csv"),
#           row.names = FALSE)

cvt <- as.data.table(cvt)
rm(cvt_data)

#' 
## ------------------------------------------------------------------------
cvt[, Conc_SD := Value_SD]
cvt[, Conc_SD_Dose := Value_SD_Dose]

#' 
#' 
#' # Evaluate fits
#' 
#' ## Table of all fit options
#' 
#' This table defines all combinations of fitting options that were tried.
#' 
## ------------------------------------------------------------------------
fit_opts <- 
  as.data.table(
    expand.grid(
  rescale_time = c(FALSE, TRUE),
  fit_conc_dose = c(FALSE, TRUE),
  fit_log_conc = c(FALSE, TRUE))
  )


#' 
#' # Evaluate residuals for each set of fitting options
#' 
## ------------------------------------------------------------------------
eval_DT <- fit_opts[,
                  {
                    print(
                      paste("fit_log_conc = ",
                        fit_log_conc,
                        "; fit_conc_dose = ",
                        fit_conc_dose,
                        "; rescale_time = ",
                        rescale_time)
                    )
                    
  #get fitting outputs file names
  fit_outputs <- dir(path = "../inst/ext",
      pattern = paste0("PK_fit_table_",
                       ".+",
        "_fit_log_conc_",
                        fit_log_conc, "_",
                        "_fit_conc_dose_",
                        fit_conc_dose, "_",
                        "_rescale_time_",
                        rescale_time, "_"
                       ))
  
  fit_outputs <- paste0("../inst/ext/",
                        fit_outputs)
  
  cat(paste0(paste(fit_outputs,
              collapse = "\n"), "\n"))
  
  #read in fitting outputs
  #flat model
  fit_flat <- fread(grep(x = fit_outputs,
                         pattern = "flat",
                         value = TRUE))
  #1compartment model
  fit_1comp <- fread(grep(x = fit_outputs,
                         pattern = "1compartment",
                         value = TRUE))
  #2compartment model
   fit_2comp <- fread(grep(x = fit_outputs,
                         pattern = "2compartment",
                         value = TRUE))
   
   fit_flat[, c("fit_conc_dose",
                "fit_log_conc",
                "rescale_time") := .(fit_conc_dose,
                                     fit_log_conc,
                                     rescale_time)]
   
    fit_1comp[, c("fit_conc_dose",
                "fit_log_conc",
                "rescale_time") := .(fit_conc_dose,
                                     fit_log_conc,
                                     rescale_time)]
    
       fit_2comp[, c("fit_conc_dose",
                "fit_log_conc",
                "rescale_time") := .(fit_conc_dose,
                                     fit_log_conc,
                                     rescale_time)]
   
   cat("Evaluating residuals...\n")
   
   #evaluate residuals
   this_eval_DT <- evaluate_resids(cvt_pre = cvt,
                         fit_flat = fit_flat,
                         fit_1comp = fit_1comp,
                         fit_2comp = fit_2comp,
                         fit_conc_dose = fit_conc_dose,
                         fit_log_conc = fit_log_conc,
                         optimx_args = list(
                           "method" = "bobyqa",
                           "itnmax" = 1e6,
                           "control" = list("kkt" = FALSE)
                         ),
                         verbose = TRUE)
   
     write.csv(this_eval_DT,
                  paste("../inst/ext/evaluate_residuals_",
                        "_fit_log_conc_",
                        fit_log_conc, "_",
                        "fit_conc_dose_",
                        fit_conc_dose, "_",
                        "rescale_time_",
                        rescale_time, "_",
                        timestamp,
                        ".csv",
                        sep=""), row.names = FALSE)
     
     this_eval_DT
                  },
  by = c("fit_conc_dose", "fit_log_conc", "rescale_time")]

write.csv(eval_DT,
          paste0("../inst/ext/eval_resid_",
          timestamp,
          ".csv"),
          row.names = FALSE)

#' 
#' Evalaute TK stats
#' 
## ------------------------------------------------------------------------
TKstats_DT <- fit_opts[,
                  {
                    print(
                      paste("fit_log_conc = ",
                        fit_log_conc,
                        "; fit_conc_dose = ",
                        fit_conc_dose,
                        "; rescale_time = ",
                        rescale_time)
                    )
                    
  #get fitting outputs file names
  fit_outputs <- dir(path = "../inst/ext",
      pattern = paste0("PK_fit_table_",
                       ".+",
        "_fit_log_conc_",
                        fit_log_conc, "_",
                        "_fit_conc_dose_",
                        fit_conc_dose, "_",
                        "_rescale_time_",
                        rescale_time, "_"
                       ))
  
  fit_outputs <- paste0("../inst/ext/",
                        fit_outputs)
  
  cat(paste0(paste(fit_outputs,
              collapse = "\n"), "\n"))
  
  #read in fitting outputs
  #flat model
  fit_flat <- fread(grep(x = fit_outputs,
                         pattern = "flat",
                         value = TRUE))
  #1compartment model
  fit_1comp <- fread(grep(x = fit_outputs,
                         pattern = "1compartment",
                         value = TRUE))
  #2compartment model
   fit_2comp <- fread(grep(x = fit_outputs,
                         pattern = "2compartment",
                         value = TRUE))
   
   fit_flat[, c("fit_conc_dose",
                "fit_log_conc",
                "rescale_time") := .(fit_conc_dose,
                                     fit_log_conc,
                                     rescale_time)]
   
    fit_1comp[, c("fit_conc_dose",
                "fit_log_conc",
                "rescale_time") := .(fit_conc_dose,
                                     fit_log_conc,
                                     rescale_time)]
    
       fit_2comp[, c("fit_conc_dose",
                "fit_log_conc",
                "rescale_time") := .(fit_conc_dose,
                                     fit_log_conc,
                                     rescale_time)]
       
       cat("Evaluating TK statistics...\n")
        #evaluate TK statistics & merge with postprocessed fits
   this_TKstats_DT <- evaluate_tkstats(cvt_pre = cvt,
                         fit_flat = fit_flat,
                         fit_1comp = fit_1comp,
                         fit_2comp = fit_2comp,
                         group_by = c("DTXSID",
                                          "Species",
                                          "Route",
                                          "Media"),
                         dose_norm = TRUE)
   
    write.csv(this_TKstats_DT,
                  paste("../inst/ext/evaluate_tkstats_",
                        "_fit_log_conc_",
                        fit_log_conc, "_",
                        "fit_conc_dose_",
                        fit_conc_dose, "_",
                        "rescale_time_",
                        rescale_time, "_",
                        timestamp,
                        ".csv",
                        sep=""), row.names = FALSE)
     
     this_TKstats_DT
                  },
  by = c("fit_conc_dose", "fit_log_conc", "rescale_time")]

write.csv(TKstats_DT,
          paste0("../inst/ext/eval_TKstats_",
          timestamp,
          ".csv"),
          row.names = FALSE)
   

#' 
#' 
