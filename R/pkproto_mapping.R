#' New mapping
#' @param mapping A [ggplot2::aes()] call that maps variable names in the
#'  original data to the harmonized `invivoPKfit` variable names.
#' @param ... Additional arguments. Currently unused.
#' @return An object of class `uneval` containing the mapping -- see
#'  [ggplot2::aes()] for details.
#' @export
#' @author Caroline Ring
mapping <- function(
    mapping = ggplot2::aes(
      Chemical = analyzed_chem_dtxsid,
      Chemical_Name = analyzed_chem_name_original,
      DTXSID = analyzed_chem_dtxsid,
      CASRN = analyzed_chem_casrn,
      Species = species,
      Reference = fk_extraction_document_id,
      Media = conc_medium_normalized,
      Route = administration_route_normalized,
      Dose = invivPK_dose_level,
      Dose.Units = "mg/kg",
      Subject_ID = fk_subject_id,
      Series_ID = fk_series_id,
      Study_ID = fk_study_id,
      ConcTime_ID = conc_time_id,
      N_Subjects = n_subjects_normalized,
      Weight = weight_kg,
      Weight.Units = "kg",
      Time = time_hr,
      Time.Units = "hours",
      Value = invivPK_conc,
      Value.Units = "mg/L",
      Value_SD = invivPK_conc_sd,
      LOQ = invivPK_loq
    ),
...
) {
  class(mapping) <- c(class(mapping), # class "uneval" by default
                      "pkproto")
return(mapping)
}
