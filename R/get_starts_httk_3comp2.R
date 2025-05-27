#' Get starting values for 1-compartment model with specific clearance
#'
#' Derive starting values for 1-compartment model parameters from available data
#'
#' This function is called internally by [get_params_1comp_cl()] and should
#' generally not be called directly by the user.
#'
#' The full set of model parameters for the 1-compartment model includes `Vdist`,
#' `kelim`, `kgutabs`, `Fgutabs`, and `Rblood2plasma`.
#' However, in this version of the one-compartment model,
#' we use the liver, alveolar glomerular flow rates (`Q_totli`, `Q_alv`, and `Q_gfr`,
#' respectively) provided by [httk].
#' Other parameters provided by [httk] include:
#' `Fup`, `Clint`, `Kblood2air`, `Rblood2plasma`, `kgutabs`, and `Fgutabs`.
#'
#' `Vdist` is calculated from estimated total clearance and `kelim`, which is
#' calculated from the data.
#'
#' The numerical optimizer requires starting guesses for the value of each
#' parameter to be estimated from the data.
#'
#' These are intended to be *very* rough starting guesses, so the algorithm here
#' is extremely naive. This function is not itself intended to produce valid
#' estimates for any of the model parameters, and it is highly unlikely to do so.
#'
#' @section Estimation of `kelim`:
#'
#' First, data are filtered to exclude any non-detects.
#'
#' Then, data are split by route of administration, into an IV data set and an oral data
#' set. (It is possible that either IV or oral data may not be
#' available for a chemical.)
#' If IV data exist, then only IV data are used to derive starting estimates for
#' `kelim`, even if oral data also exist.
#'
#' If only oral data exist, then the oral data are used to derive a starting
#' estimate for `kelim`.
#'
#' Whichever data set is used (IV or oral), the starting value for `kelim` is
#' derived by assuming that the range of observed time values in the data set
#' spans two elimination half-lives. This implies that the elimination half-life
#' is equal to the midpoint of observed time values, and that the starting value
#' for the elimination time constant `kelim` is therefore `log(2)` divided by the
#' midpoint of observed time values.
#'
#' Of course, this assumption is unlikely to be correct. However, we hope that it
#' will yield a starting guess for `kelim` that is at least on the right order of
#' magnitude.
#'
#' @section Starting value for `Vdist`:
#'
#' Using a calculated value for total clearance, `Cltot`, `Vdist` is estimated
#' by dividing this by the estimation of `kelim`.
#'
#' @inheritParams get_starts_flat
#' @param this_chemical A character vector naming the chemical for calculations in `httk`.
#' @param this_species A character vector naming the species for calculations in `httk`.
#' @param restrictive A boolean value determinining whether to assume restrictive
#'   or non-restrictive clearance when getting starting values.
#'
#' @return The same `data.frame` as `par_DF`, with an additional variable
#'  `starts` containing the derived starting value for each parameter. If a
#'  parameter cannot be estimated from the available data, then its starting value
#'  will be `NA_real_`
#' @import httk
#' @author Caroline Ring
#' @family httk model functions
#' @family get_starts functions
#' @family built-in model functions
#'
get_starts_httk_3comp2 <- function(data,
                                   par_DF,
                                   this_chemical,
                                   this_species,
                                   restrictive) {

  parm_3comp2 <- tryCatch(
    expr = {
      httk::parameterize_3comp2(
        dtxsid = unique(this_chemical),
        species = unique(this_species),
        default.to.human = TRUE,
        restrictive.clearance = restrictive
      ) |>
        suppressWarnings() |>
        suppressMessages()
    }, error = function(e) {
      message("Error: ", e)
      if (interactive()) {
        response <- readline(
          prompt = paste0(
            "There has been an error, ",
            "substitute with starting parameters for Bisphenol A?"
          )
        ) |>
          tolower() |>
          trimws()
        if (startsWith(response, 'y')) {
          httk::parameterize_3comp2(
            dtxsid = "DTXSID7020182",
            species = this_species,
            default.to.human = TRUE,
            restrictive.clearance = restrictive
          ) |>
            suppressWarnings() |>
            suppressMessages()
        } else {
          # Early return with all values set to NA_real_
          c(
            "BW" = NA_real_,
            "Caco2.Pab" = NA_real_,
            "Caco2.Pab.dist" = NA_real_,
            "Clint" = 0,
            "Clint.dist" = NA_real_,
            "Clmetabolismc" = NA_real_,
            "Funbound.plasma" = 0,
            "Funbound.plasma.dist" = NA_real_,
            "Funbound.plasma.adjustment" = NA_real_,
            "Fabsgut" = NA_real_,
            "Fhep.assay.correction" = NA_real_,
            "hematocrit" = NA_real_,
            "Kgut2pu" = NA_real_,
            "Krbc2pu" = NA_real_,
            "kgutabs" = NA_real_,
            "Kliver2pu" = NA_real_,
            "Krest2pu" = NA_real_,
            "Kblood2air" = NA_real_,
            "liver.density" = NA_real_,
            "logHenry" = NA_real_,
            "million.cells.per.gliver" = NA_real_,
            "MW" = NA_real_,
            "Pow" = NA_real_,
            "pKa_Donor" = NA_real_,
            "pKa_Accept" = NA_real_,
            "MA" = NA_real_,
            "Qcardiacc" = NA_real_,
            "Qgfrc" = NA_real_,
            "Qgutf" = NA_real_,
            "Qliverf" = NA_real_,
            "Qalvc" = NA_real_,
            "Rblood2plasma" = NA_real_,
            "Vgutc" = NA_real_,
            "Vliverc" = NA_real_,
            "Vrestc" = NA_real_
          )
        }
      }
    }
  )

  starts <- parm_3comp2

  if (all(sapply(parm_3comp2, is.na))) {

    starts["pKa_Donor"] = c(pKa_Donor = ' ')
    starts["pKa_Accept"] = c(pKa_Accept = ' ')

    par_DF$start <- as.numeric(starts[par_DF$param_name])
    return(par_DF)
  }

  units_pKa <- "logarithmic"
  # Check the value of pKa_Donor and pKa_Accept
  tmp_pKa_Donor <- string2num(starts[["pKa_Donor"]])
  len_pKa_Donor <- length(tmp_pKa_Donor)
  nom_pKa_Donor <- paste0("pKa_Donor", "_", seq_len(len_pKa_Donor))

  df_pKa_Donor <- data.frame("param_name" = nom_pKa_Donor,
                          "param_units" = units_pKa,
                          "optimize_param" = FALSE,
                          "use_param" = TRUE,
                          "lower_bound" = NA_real_,
                          "upper_bound" = NA_real_,
                          "start" = tmp_pKa_Donor)

  tmp_pKa_Accept <- string2num(starts[["pKa_Accept"]])
  len_pKa_Accept <- length(tmp_pKa_Accept)
  nom_pKa_Accept <- paste0("pKa_Accept", "_", seq_len(len_pKa_Accept))

  df_pKa_Accept <- data.frame("param_name" = nom_pKa_Accept,
                          "param_units" = units_pKa,
                          "optimize_param" = FALSE,
                          "use_param" = TRUE,
                          "lower_bound" = NA_real_,
                          "upper_bound" = NA_real_,
                          "start" = tmp_pKa_Accept)
  # Note: Even when the pKa = ' ' = NA, there will be a pKa_[Donor/Accept]_1

  # Assemble the entire parameter data.frame
  par_DF$start <- as.numeric(starts[par_DF$param_name])
  par_DF <- rbind(par_DF, df_pKa_Donor, df_pKa_Accept)
  rownames(par_DF) <- NULL

  return(par_DF)
}
