#' New mapping
#' @param mapping A [ggplot2::aes()] call that maps variable names in the original data to the harmonized `invivoPKfit` variable names.
#' @param ... Additional arguments. Currently unused.
#' @return An object of class `uneval` containing the mapping -- see [ggplot2::aes()] for details.
#' @export
#' @author Caroline Ring
mapping <- function(
    mapping = ggplot2::aes(
      Chemical = analyte_dtxsid,
      Chemical_Name = analyte_name_original,
      DTXSID = analyte_dtxsid,
      CASRN = analyte_casrn,
      Species = species,
      Reference = document_id,
      Media = conc_medium_normalized,
      Route = administration_route_normalized,
      Dose = dose_level_normalized,
      Dose.Units = "mg/kg",
      Subject_ID = subject_id,
      Series_ID = series_id,
      Study_ID = study_id,
      ConcTime_ID = conc_time_id,
      N_Subjects = n_subjects_in_series,
      Weight = weight_kg,
      Weight.Units = "kg",
      Time = time_hr,
      Time.Units = "hours",
      Value = conc,
      Value.Units = "mg/L",
      Value_SD = conc_sd_normalized,
      LOQ = loq
    ),
...
){
  class(mapping) <- c(class(mapping), #class "uneval" by default
                      "pkproto")
return(mapping)
}
