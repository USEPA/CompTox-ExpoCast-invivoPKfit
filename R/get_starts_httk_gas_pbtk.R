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
#' @author Gilberto Padilla Mercado
#' @family httk model functions
#' @family get_starts functions
#' @family built-in model functions
#'
get_starts_httk_gas_pbtk <- function(data,
                                     par_DF,
                                     this_chemical,
                                     this_species,
                                     restrictive) {
  this_chemical <- unique(this_chemical)
  this_species <- unique(this_species)
  stopifnot(all(lengths(list(this_chemical, this_species)) == 1))

  parm_gas_pbtk <- tryCatch(
    expr = {
      httk::parameterize_gas_pbtk(
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
        if (startsWith(response, "y")) {
          httk::parameterize_gas_pbtk(
            dtxsid = "DTXSID7020182",
            species = this_species,
            default.to.human = TRUE,
            restrictive.clearance = restrictive
          ) |>
            suppressWarnings() |>
            suppressMessages()
        } else {
          # Early return with all values set to NA_real_
          list(
            "BW" = NA_real_,
            "Caco2.Pab" = NA_real_,
            "Caco2.Pab.dist" = NA_character_,
            "Clint" = NA_real_,
            "Clint.dist" = NA_character_,
            "Clmetabolismc" = NA_real_,
            "Funbound.plasma" = NA_real_,
            "Funbound.plasma.dist" = NA_real_,
            "Funbound.plasma.adjustment" = NA_real_,
            "Fabsgut" = NA_real_,
            "Fhep.assay.correction" = NA_real_,
            "hematocrit" = NA_real_,
            "Kgut2pu" = NA_real_,
            "Krbc2pu" = NA_real_,
            "kgutabs" = NA_real_,
            "Kkidney2pu" = NA_real_,
            "Klung2pu" = NA_real_,
            "km" = NA_real_,
            "Kmuc2air" = NA_real_,
            "Kliver2pu" = NA_real_,
            "Krest2pu" = NA_real_,
            "Kblood2air" = NA_real_,
            "kUrtc" = NA_real_,
            "liver.density" = 1.05,
            "logHenry" = NA_real_,
            "million.cells.per.gliver" = 110,
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
            "Qkidneyf" = NA_real_,
            "Qlungf" = NA_real_,
            "Rblood2plasma" = NA_real_,
            "Vgutc" = NA_real_,
            "Vliverc" = NA_real_,
            "Vartc" = NA_real_,
            "Vkidneyc" = NA_real_,
            "Vlungc" = NA_real_,
            "vmax" = 0,
            "Vmucc" = NA_real_,
            "Vvenc" = NA_real_,
            "Vrestc" = NA_real_,
            "KFsummary" = NA_real_,
            "Fprotein.plasma" = NA_real_,
            "fabs.oral" = NA_real_,
            "Qgut_" = NA_real_,
            "Qintesttransport" = NA_real_
          )
        }
      }
    }
  )

  starts <- parm_gas_pbtk
  # Following are required for calculations of 'extra' parameters.
  stopifnot(all(!is.na(starts[c("Caco2.Pab", "Funbound.plasma", "krbc2pu", "BW", "Qgutf", "Qcardiacc")])))

  fabs.oral <- httk::calc_fabs.oral(list(Caco2.Pab = starts[["Caco2.Pab"]]),
                                    species = this_species)
  Fprotein.plasma <- httk::physiology.data[
    which(httk::physiology.data[, "Parameter"] == "Plasma Protein Volume Fraction"),
    which(tolower(colnames(httk::physiology.data)) == tolower(this_species))
  ]
  Kint <- 1 - Fprotein.plasma +
    (0.37 / starts[["Funbound.plasma"]] - (1 - Fprotein.plasma))
  KFsummary <- starts[["Krbc2pu"]] / Kint
  Qintesttransport <- 0.1 * (starts[["BW"]] / 70)^(3/4)

  peff <- 10^(0.4926 * log10(starts[["Caco2.Pab"]]) - 0.1454)
  Asi <- 0.66 * starts[["BW"]] / 70
  if (this_species == "rat") {
    peff <- max(0, (peff + 0.1815)/(1.039 * 10))
    Asi <- 71/(100^2)
  }
  CLperm <- peff * Asi * 36000
  Qvilli <- (18 / (38.7/1000 * 60 * 15.8757)) *
    starts[["Qcardiacc"]] * starts[["Qgutf"]] * starts[["BW"]]^(3/4)

  Qgut_ <- Qvilli * CLperm / (Qvilli + CLperm)


  starts[["KFsummary"]] <- signif(KFsummary, 8)
  starts[["Fprotein.plasma"]] <- signif(Fprotein.plasma, 8)
  starts[["fabs.oral"]] <- signif(fabs.oral, 8)
  starts[["Qgut_"]] <- signif(Qgut_, 8)
  starts[["Qintesttransport"]] <- signif(Qintesttransport, 8)

  # To pass check, values must not be NA, resolve in downstream functions
  if (is.na(starts[["Caco2.Pab.dist"]])) {
    starts[["Caco2.Pab.dist"]] <- ""
  }
  if (is.na(starts[["Clint.dist"]])) {
    starts[["Clint.dist"]] <- ""
  }
  if (is.na(starts[["Funbound.plasma.dist"]])) {
    starts[["Funbound.plasma.dist"]] <- ""
  }

  # Assemble the entire parameter data.frame
  par_DF$start <- starts[par_DF$param_name]
  rownames(par_DF) <- NULL

  return(par_DF)
}
