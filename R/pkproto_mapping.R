#' New mapping
#' @param mapping A [ggplot2::aes()] call that maps variable names in the original data to the harmonized `invivoPKfit` variable names.
#' @param ... Additional arguments. Currently unused.
#' @return An object of class `uneval` containing the mapping -- see [ggplot2::aes()] for details.
#' @export
#' @author Caroline Ring
mapping <- function(
    mapping = ggplot2::aes(
  Chemical = chemicals_analyzed.dsstox_substance_id,
                                 Chemical_Name = series.analyte_name_original,
                                 DTXSID = chemicals_analyzed.dsstox_substance_id,
                                 CASRN = chemicals_analyzed.dsstox_casrn,
                                 Species = subjects.species_harmonized,
                                 Reference = as.character(
                                   ifelse(
                                     is.na(
                                       documents_reference.id
                                     ),
                                     documents_extraction.id,
                                     documents_reference.id
                                   )
                                 ),
                                 Media = series.conc_medium_normalized,
                                 Route = studies.administration_route_normalized,
                                 Dose = studies.dose_level_normalized_corrected,
                                 Dose.Units = "mg/kg",
                                 Subject = subjects.id,
                                 Series_ID = series.id,
                                 Study_ID = studies.id,
                                 ConcTime_ID = conc_time_values.id,
                                 N_Subjects =  series.n_subjects_in_series,
                                 Weight = subjects.weight_kg,
                                 Weight.Units = "kg",
                                 Time = conc_time_values.time_hr,
                                 Time.Units = "hours",
                                 Value = conc_time_values.conc,
                                 Value.Units = "mg/L",
                                 LOQ = series.loq_normalized,
                                 Value_SD  = conc_time_values.conc_sd_normalized
),
...
){
  class(mapping) <- c(class(mapping), #class "uneval" by default
                      "pkproto")
return(mapping)
}
