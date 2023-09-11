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

#' CvTdb data
#'
#' Concentration vs. time data from CvTdb
#'
#' This is concentration vs. time data from CvTdb, most recently downloaded as
#' of the date in [`cvt_download_date`].
#'
#' These data have been filtered to retain only oral and intravenous
#' administration, and only measurements in blood and plasma. They have also
#' been filtered to retain only observations where the same chemical was both
#' administered and measured in blood/plasma (i.e., excluding observations where
#' a metabolite was measured).
#'
#' A `data.frame` with 14805 rows and 105 variables.
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
