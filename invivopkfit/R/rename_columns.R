#'Harmonize variable names.
#'
#'Renames concentration vs. time data to use harmonized variable names.
#'
#'# Default values for `names_list`
#'
#' ```
#' list(
#'"Compound_Dosed" = "studies.test_substance_name_original",
#'"DTXSID_Dosed" = "chemicals_dosed.dsstox_substance_id",
#'"CAS_Dosed" = "chemicals_dosed.dsstox_casrn",
#'"Compound_Analyzed" = "series.analyte_name_original",
#' "DTXSID_Analyzed" = "chemicals_analyzed.dsstox_substance_id",
#' "CAS_Analyzed" = "chemicals_analyzed.dsstox_casrn",
#' "Reference" = "documents_reference.id",
#' "Extraction" = "documents_extraction.id",
#' "Species" = "subjects.species",
#' "Weight" ="subjects.weight_kg",
#' "Weight.Units" = NULL,
#' "Dose" = "studies.dose_level_normalized",
#' "Dose.Units" = NULL,
#' "Time" = "conc_time_values.time_hr",
#' "Time.Units" = NULL,
#' "Media" = "series.conc_medium_normalized",
#' "Value" = "conc_time_values.conc",
#' "Value.Units" = NULL,
#' "Route" = "studies.administration_route_normalized",
#' "LOQ" = "series.loq",
#' "Subject" = "subjects.id",
#' "N_Subjects" = "series.n_subjects_in_series",
#' "Value_SD" = "conc_time_values.conc_sd")
#' ```
#'
#'# Default values for `defaults_list`
#'
#' ```
#' list(
#' "Weight.Units" = "kg",
#' "Dose.Units" = "mg/kg",
#' "Time.Units" = "h",
#' "Value.Units" = "mg/L")
#' ```
#'
#'@param data.set A \code{data.frame} of concentration-time data. Must contain
#'  the variables specified in `names_list`.
#'@param names_list A named list of format `list(new_name = old_name)`, where
#'  the names are the new variable names, and the list elements (character) are
#'  the old variable names (in `data.set`).  Each old variable name will be
#'  renamed to its corresponding new variable name. If a named list element is
#'  NULL, then there is no old variable corresponding to that new variable. In
#'  this case, a new variable will be added, with default value as given in
#'  `defaults_list`; if no default is given in `defaults_list`, the new variable
#'  will be filled with `NA_character`. See Details for default value.
#'@param defaults_list A named list where the names are the new variable names,
#'  and the values are the default values to fill for each new variable, if the
#'  corresponding old variable name is NULL in `names_list`. See Details for
#'  default value.
#'
#'@return A `data.frame` that contains the variables specfiied in `names_list`,
#'  but renamed to use the harmonized names given as the names of `names_list`.
#'@author John Wambaugh, Chris Cook, Caroline Ring

rename_columns <- function(data.set,
                           names_list =list(
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
                             "Dose" = "studies.dose_level_normalized",
                             "Dose.Units" = NULL,
                             "Time" = "conc_time_values.time_hr",
                             "Time.Units" = NULL,
                             "Media" = "series.conc_medium_normalized",
                             "Value" = "conc_time_values.conc",
                             "Value.Units" = NULL,
                             "Route" = "studies.administration_route_normalized",
                             "LOQ" = "series.loq",
                             "Subject" = "subjects.id",
                             "N_Subjects" = "series.n_subjects_in_series",
                             "Value_SD" = "conc_time_values.conc_sd"),
                           defaults_list =   list(
                             "Weight.Units" = "kg",
                             "Dose.Units" = "mg/kg",
                             "Time.Units" = "h",
                             "Value.Units" = "mg/L")) {

  data.set <- as.data.frame(data.set)

  #if any old names are not present in data.set,
  #throw an errror
  missing_items <- names_list[!(names_list %in% names(data.set)) &
                                sapply(names_list, function(x) !is.null(x))]
  if(length(missing_items)>0){
    warning(paste0("invivopkfit::rename_columns():\n",
                  "The following items in names_list are columns not present in data.set:\n",
                 paste(
                   paste(names(missing_items),
                       missing_items,
                       sep = " = "),
                       collapse = "\n"
                   ),
                 "\nThey will be treated as though they were NULL."))
    names_list[!names_list %in% names(data.set)] <- list(NULL)
  }

  #new column names
  new_names <- names(names_list[sapply(names_list, function(x) !is.null(x))])
  #old column names
  old_names <- unlist(names_list[sapply(names_list, function(x) !is.null(x))])
  #columns to add: had NULL instead of an old_name
  add_names <- names(names_list)[sapply(names_list, is.null)]

  #pull the specified old names
  dat <- data.set[old_names]

  #rename to new names
  names(dat) <- new_names

  #add the columns that had no old names
  dat[add_names] <- as.list(rep(NA_character_, length(add_names)))

  #fill with default values if provided
  dat[names(defaults_list)] <- defaults_list

  return(dat)
}
