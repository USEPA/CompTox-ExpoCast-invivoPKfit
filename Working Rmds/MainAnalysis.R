## ----setup, include = FALSE-----------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(results = TRUE)
knitr::opts_chunk$set(message = FALSE)


## ---- include = FALSE-----------------------------------------------------
library(tidyverse)
library(data.table)
library(readxl)
devtools::load_all("../invivopkfit")


## ---- load_cvt, eval = TRUE-----------------------------------------------
cvt_data <- fread("../inst/ext/cvt_data_2023_01_05.csv")


## -------------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y_%m_%d")


## -------------------------------------------------------------------------
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


## -------------------------------------------------------------------------
defaults_list <- list(
                      "Weight.Units" = "kg",
                      "Dose.Units" = "mg/kg",
                      "Time.Units" = "hours",
                      "Value.Units" = "mg/L")


## -------------------------------------------------------------------------
cvt <- preprocess_data(data.set = cvt_data,
                       names_list =names_list,
                       defaults_list =   defaults_list,
                       ratio_conc_to_dose = 1,
                       calc_loq_factor = 0.45,
                       routes_keep = c("po", "iv"),
                       media_keep = c("blood", "plasma"),
                       impute_loq = TRUE,
                       impute_sd = TRUE,
                       study_def = c("DTXSID", "Species", "Reference"),
                       suppress.messages = TRUE)

#save pre-processed data
write.csv(cvt,
          paste0("../inst/ext/cvt_preprocessed_",
                 timestamp,
          ".csv"),
          row.names = FALSE)


## -------------------------------------------------------------------------
for(this_method in c("bobyqa", "L-BFGS-B")){
  for(this_rescale_time in c(FALSE, TRUE)){
    for(this_fit_conc_dose in c(FALSE, TRUE)){
      for(this_fit_log_conc in c(FALSE, TRUE)){
          for (this_model in c("flat",
                               "1compartment",
                               "2compartment")){
            cat(paste0("method = ", this_method, "\n",
                         "fit_log_conc = ", this_fit_log_conc, "\n",
                         "fit_conc_dose = ", this_fit_conc_dose, "\n",
                         "rescale_time = ", this_rescale_time, "\n",
                                     "model = ", this_model, "\n"))
            system.time(
              PK.fit.table <- fit_all(
                data.set = cvt,
                model = this_model,
                modelfun = "analytic",
                preprocess = FALSE,
                fit_conc_dose = this_fit_conc_dose,
                fit_log_conc = this_fit_log_conc,
                rescale_time = this_rescale_time,
                optimx_args = list(
                  method = this_method,
                  itnmax = 1e6,
                  control = list(kkt = FALSE,
                                 factr = 1e7)
                ),
                suppress.messages = FALSE
              )
            )
            
            write.csv(PK.fit.table,
                      paste0("../inst/ext/PK_fit_table_",
                            this_model, "_",
                            "fit_log_conc_",
                            this_fit_log_conc, "_",
                            "fit_conc_dose_",
                            this_fit_conc_dose, "_",
                            "rescale_time_",
                            this_rescale_time, "_",
                            "method_", this_method, "_",
                            timestamp,
                            ".csv"), row.names = FALSE)
            
          } #end loop over models
      } #end loop over this_fit_log_conc
    } #end loop over this_fit_conc_dose
  } #end loop over this_rescale_time
}#end loop over this_method

