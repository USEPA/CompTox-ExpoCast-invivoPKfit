library(data.table)
devtools::load_all("invivopkfit")
library(dplyr)


cvt <- fread("inst/ext/cvt_data_2023_01_05.csv",
             na.strings = c("", "NA"))

cvt <- cvt[chemicals_analyzed.dsstox_substance_id == chemicals_dosed.dsstox_substance_id &
           studies.administration_route_normalized %in% c("oral", "iv") &
           series.conc_medium_normalized %in% c("blood", "plasma") &
             (studies.dose_level_normalized > 0) %in% TRUE &
             (conc_time_values.time_hr >= 0) %in% TRUE, ]

cvt_sub <- cvt[chemicals_analyzed.dsstox_substance_id %in% "DTXSID3061635"]


my_pk <- pk(data = cvt_sub,
            mapping = ggplot2::aes(Chemical = chemicals_analyzed.dsstox_substance_id,
                                   DTXSID = chemicals_analyzed.dsstox_substance_id,
                                   Chemical_Name = chemicals_analyzed.preferred_name,
                                   CASRN = chemicals_analyzed.dsstox_casrn,
                                   Species = subjects.species,
                                   Reference = as.character(ifelse(is.na(documents_reference.id),
                                                      documents_extraction.id,
                                                      documents_reference.id)),
                                   Media = series.conc_medium_normalized,
                                   Route = studies.administration_route_normalized,
                                   Dose = studies.dose_level_normalized,
                                   Dose.Units = "mg/kg",
                                   Subject = subjects.id,
                                   N_Subjects =  series.n_subjects_in_series,
                                   Weight = subjects.weight_kg,
                                   Weight.Units = "kg",
                                   Time = conc_time_values.time_hr,
                                   Time.Units = "hours",
                                   Value = conc_time_values.conc,
                                   Value.Units = "mg/L",
                                   LOQ = series.loq_normalized,
                                   Value_SD  = conc_time_values.conc_sd_normalized,
                                   Conc = pmax(conc_time_values.conc,
                                               series.loq_normalized,
                                               na.rm = TRUE),
                                   Detect = ifelse(is.na(conc_time_values.conc) |
                                                     conc_time_values.conc <= series.loq_normalized,
                                                   FALSE,
                                                   TRUE),
                                   Conc_SD = conc_time_values.conc_sd_normalized)
                                   )

my_pk <- preprocess_data(my_pk)
