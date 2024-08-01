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
                                          document_id = "Reference",
                                          conc_medium_normalized = "Media",
                                          administration_route_normalized = "Route",
                                          dose_level_corrected = "Dose",
                                          subject_id = "Subject_ID",
                                          series_id = "Series_ID",
                                          study_id = "Study_ID",
                                          conc_time_id = "ConcTime_ID",
                                          n_subjects_in_series = "N_Subjects",
                                          weight_kg = "Weight",
                                          time_hr = "Time",
                                          conc = "Value",
                                          conc_sd_normalized = "Value_SD",
                                          loq = "LOQ")) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame!")
  }

  new_data <- dplyr::rename(data, tidyr::any_of(cvt_LUT))

  return(new_data)
}
