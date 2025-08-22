#' Allowable time units
#'
#' A `character` vector of allowable units for time variables.
#'
#' These are the time units understood by [lubridate::period()] and [lubridate::duration()].
"time_units"

#' Time conversion table
#'
#' A `data.frame` that has the converted units from "time_units"
#'
"time_conversions"

#' 1-compartment model
#'
#' The `pk_model` object defining the 1-compartment model.
#'
#' A `pk_model` object: under the hood, a `list` object with named elements
#' corresponding to the arguments of [pk_model()]. See that function
#' documentation for the definition of each element.
#'
#' See [cp_1comp()] for the function that predicts blood/plasma concentration for a bolus dose (oral or IV).
#'
#' See [auc_1comp()] for the function that predicts area under the concentration-time curve.
#'
#' See [tkstats_1comp()] for the function that calculates summary toxicokinetic statistics from 1-compartment model parameters.
#'
#' See [params_1comp()] for the function that determines bounds and starting guesses for model parameters, based on the data.
#'
"model_1comp"

#' 2-compartment model
#'
#' The `pk_model` object defining the 2-compartment model.
#'
#' A `pk_model` object: under the hood, a `list` object with elements
#' corresponding to the arguments of [pk_model()]. See that function
#' documentation for the definition of each element.
#'
#' See [cp_2comp()] for the function that predicts blood/plasma concentration for a bolus dose (oral or IV).
#'
#' See [auc_2comp()] for the function that predicts area under the concentration-time curve.
#'
#' See [tkstats_2comp()] for the function that calculates summary toxicokinetic statistics from 1-compartment model parameters.
#'
#' See [params_2comp()] for the function that determines bounds and starting guesses for model parameters, based on the data.
"model_2comp"

#' Flat model
#'
#' The `pk_model` object defining the flat model.
#'
#' A `pk_model` object: under the hood, a `list` object with elements
#' corresponding to the arguments of [pk_model()]. See that function
#' documentation for the definition of each element.
#'
#' See [cp_flat()] for the function that predicts blood/plasma concentration for a bolus dose (oral or IV).
#'
#' See [auc_flat()] for the function that predicts area under the concentration-time curve.
#'
#' See [tkstats_flat()] for the function that calculates summary toxicokinetic statistics from 1-compartment model parameters.
#'
#' See [params_flat()] for the function that determines bounds and starting guesses for model parameters, based on the data.
"model_flat"

#' Gas pbtk `httk` model
#'
#' The `pk_model` object defining the "gas_pbtk" model from `httk`, with
#' recalculations of other "constant" parameters that depend on optimized parameters.
#'
#' A `pk_model` object: under the hood, a `list` object with named elements
#' corresponding to the arguments of [pk_model()]. The default set of parameters
#' to optimize are `Clint` and `Funbound.plasma`, but users can also fit partition coefficients.
#' Note that fitting may present instability during ODE solving step.
#' See that function documentation for the definition of each element.
#'
#' See [cp_httk_gas_pbtk()] for the function that predicts blood/plasma concentration for a bolus dose (oral or IV).
#'
#' See [get_params_httk_gas_pbtk()] for the function that determines bounds and starting guesses for model parameters, based on the data.
#'
"model_httk_gas_pbtk"


#' CvTdb data
#'
#' Concentration vs. time data from CvTdb
#'
#' This is concentration vs. time data from CvTdb, most recently downloaded as
#' of the date in [`cvt_date`].
#'
#' These data have been filtered to retain only oral and intravenous
#' administration, and only measurements in blood and plasma. They have also
#' been filtered to retain only observations where the same chemical was both
#' administered and measured in blood/plasma (i.e., excluding observations where
#' a metabolite was measured).
#'
#' @format A `data.frame` with 13937 rows and 61 variables:
#' \describe{
#'   \item{conc_time_id}{Unique database identifier for each CvT observation.}
#'   \item{fk_series_id}{Unique database identifier for experimental series.}
#'   \item{time_original}{Timepoint in original units.}
#'   \item{time_hr}{Timepoint in hours.}
#'   \item{conc_original}{Concentration in original units.}
#'   \item{conc_sd_original}{Standard deviation of concentration in original units.}
#'   \item{conc}{Concentration in normalized units.}
#'   \item{conc_sd}{Standard deviation of concentration in normalized units.}
#'   \item{fk_analyzed_chemical_id}{Unique database identifier for analyte.}
#'   \item{analyzed_chem_dtxsid}{DTXSID of chemical analyte.}
#'   \item{analyzed_chem_name_original}{Original analyte name.}
#'   \item{analyzed_chem_casrn}{CASRN of chemical analyte.}
#'   \item{analyzed_chem_name}{Preferred name for analyte.}
#'   \item{time_units_original}{Original time units.}
#'   \item{conc_units_original}{Original concentration units.}
#'   \item{conc_units_normalized}{Normalized concentration units.}
#'   \item{conc_unit_norm_factor}{Ratio of conc/conc_original}
#'   \item{loq}{Level of quantification.}
#'   \item{loq_units}{Units for loq.}
#'   \item{n_subjects_in_series}{Number of subjects in each series.}
#'   \item{radiolabeled}{Answers whether this observation comes from a radiolabelling or isotope tracing experiment.}
#'   \item{fk_study_id}{Unique database identifier for each study.}
#'   \item{administration_route_normalized}{Route of exposure/administration, either oral or iv.}
#'   \item{fk_dosed_chemical_id}{Unique database identifier for dosed chemical.}
#'   \item{dosed_chem_dtxsid}{DTXSID of dosed chemical.}
#'   \item{dosed_chem_name_original}{Original dosed chemical name.}
#'   \item{dosed_chem_casrn}{CASRN of dosed chemical.}
#'   \item{dosed_chem_name}{Preferred name for dosed chemical.}
#'   \item{dose_volume}{Volume of dose.}
#'   \item{dose_volume_units}{Units for dose_volume.}
#'   \item{dose_vehicle}{If available, specifies what vehicle was used CvT experiment.}
#'   \item{dose_duration}{Duration of the dose, if available.}
#'   \item{dose_duration_units}{Units for dose_duration.}
#'   \item{dose_frequency}{Frequency of dosing, these should all be 1 for a single bolus.}
#'   \item{fasting_period}{If available, describes the fasting period for subjects.}
#'   \item{dose_level_normalized}{Dose levels in normalized units.}
#'   \item{dose_level_original}{Dose levels in original units.}
#'   \item{dose_level_units_original}{Units for dose_level_original.}
#'   \item{conc_medium_normalized}{Standardized media names, blood or plasma.}
#'   \item{conc_medium_original}{Original media names.}
#'   \item{fk_subject_id}{Unique database identifier for each subject.}
#'   \item{weight_kg}{Subject weight, in kilograms.}
#'   \item{species}{Subject species.}
#'   \item{sex}{Subject sex, if available.}
#'   \item{age}{Subject age, if available.}
#'   \item{age_units}{Units for age.}
#'   \item{age_category}{Categories for age.}
#'   \item{fk_extraction_document_id}{Unique database identifier for documents that were curated.}
#'   \item{pmid}{Document PubMed ID.}
#'   \item{year}{Year of publication.}
#'   \item{other_study_identifier}{Alternative identifier for documents, used for NTP studies.}
#'   \item{url}{Document URL address.}
#'   \item{doi}{Document DOI.}
#'   \item{extracted}{Curation level of document.}
#'   \item{curation_set_tag}{A grouping tag for specific document extraction and curation efforts.}
#'   \item{n_subjects_normalized}{Normalized subject number.}
#'   \item{invivPK_dose_level_units}{dose_level_units used in this package, mg/kg.}
#'   \item{invivPK_conc_units}{conc_units used in this package, ug/mL}
#'   \item{invivPK_conc}{Concentrations normalized to ug/mL}
#'   \item{invivPK_dose_level}{Dose normalized to mg/kg}
#'   \item{invivPK_loq}{Level of quantification in ug/mL}
#'   \item{invivPK_loq_units}{LOQ units used in this package, ug/mL.}
#'   \item{invivPK_conc_sd}{Standard deviation of concentrations, in ug/mL.}
#' }
#'
"cvt"

#' CvTdb download date
#'
#' The most recent download date of [`cvt`] data
#'
#' A character scalar giving the date in "YYYY-MM-DD" format of the download
#' date of the data in [`cvt`] from the CvTdb database.
#'
"cvt_date"

#' CvTdb data for invivoPKfit 2.0.0 release (old)
#'
#' The CvTdb data released for the manuscript Informatics for toxicokinetics (2025).
#'
#' A data.frame with similar data to [`cvt`]
#'
"cvt_2.0.0"

#' SQL query result (current)
#'
#' This is the raw SQL query result.
#'
#' A data.frame similar to [`cvt`] and [`cvt_2.0.0`], but to create those objects
#' some values are changed to normalized values for use with invivoPKfit specifically
#' and to include experiments that may not pass filtering due to data coding
#' details, but are reasonable to include for analysis.
"cvtdb_original"

#' Status ID for initialization
#'
#' An integer status that denotes a [pk()] object has been initialized.
"status_init"

#' Status ID for preprocessing
#'
#' An integer status that denotes [preprocess_data()] has been completed.
"status_preprocess"

#' Status ID for data summary info
#'
#' An integer status that denotes [data_info()] has been completed.
"status_data_info"

#' Status ID for pre-fitting
#'
#' An integer status that denotes [prefit()] has been completed.
"status_prefit"

#' Status ID for fitting
#'
#' An integer status that denotes [fit()] has been completed.
"status_fit"
