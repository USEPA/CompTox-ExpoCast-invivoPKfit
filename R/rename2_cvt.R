#' Convert invivoPKfit output table names to the CvTdb names
#'
#' @param data A data frame from on the invivoPKfit outputs.
#' @param cvt_LUT A look-up table for the name conversions, can be constomized.
#' must be a vector.

rename2_cvt <- function(data,
                              cvt_LUT = c(analyte_dtxsid = "Chemical",
                                          analyte_dtxsid = "DTXSID",
                                          analyte_name_original = "Chemica_Name",
                                          species = "Species",
                                          fk_extraction_document_id = "Reference",
                                          conc_medium_normalized = "Media",
                                          administration_route_normalized = "Route",
                                          invivPK_dose_level = "Dose",
                                          fk_subject_id = "Subject_ID",
                                          fk_series_id = "Series_ID",
                                          fk_study_id = "Study_ID",
                                          conc_time_id = "ConcTime_ID",
                                          invivPK_subjects_corrected = "N_Subjects",
                                          weight_kg = "Weight",
                                          time_hr = "Time",
                                          invivPK_conc = "Value",
                                          invivPK_conc_sd = "Value_SD",
                                          invivPK_loq = "LOQ")) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame!")
  }

  new_data <- dplyr::rename(data, tidyr::any_of(cvt_LUT))
  new_data <- dplyr::rename_with(new_data,
                                 ~ tolower(gsub("\\.|\\/", "_", .x)),
                                 tidyselect::everything())
  return(new_data)
}
