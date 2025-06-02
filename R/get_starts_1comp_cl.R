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
#' @param restrictive A boolean value determinining whether to assume restrictive
#'   or non-restrictive clearance when getting starting values.
#'
#' @return The same `data.frame` as `par_DF`, with an additional variable
#'  `starts` containing the derived starting value for each parameter. If a
#'  parameter cannot be estimated from the available data, then its starting value
#'  will be `NA_real_`
#' @import httk
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family 1-compartment model functions
#' @family get_starts functions
#' @family built-in model functions
#'
get_starts_1comp_cl <- function(data,
                                par_DF,
                                restrictive) {
  # initialize starting values for each parameter.
  # if no IV data exist, then Vdist starting value will remain NA.
  # if no oral data exist, then Fgutabs_Vdist and Fgutabs starting values will remain NA.
  # if only one of IV or oral data exist, then Fgutabs starting value will remain NA.
  # May just make a static version of this?
  # Or combine it with another table
  kelim <- NA_real_
  kgutabs <- NA_real_
  Vdist <- NA_real_
  Fgutabs <- NA_real_
  Fup <- NA_real_
  this_species <- unique(data$Species)

  parm_gas <- tryCatch(
    expr = {
      httk::parameterize_3comp2(
        dtxsid = unique(data[["Chemical"]]),
        species = this_species,
        default.to.human = FALSE,
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
            "substitute with starting parameters for bisphenol A?"
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
          c("Q_totli" = NA_real_,
            "Q_gfr" = NA_real_,
            "Q_alv" = NA_real_,
            "Kblood2air" = NA_real_,
            "Fup" = NA_real_,
            "Clint" = NA_real_,
            "kgutabs" = NA_real_,
            "Vdist" = NA_real_,
            "Fgutabs" = NA_real_,
            "Rblood2plasma" = NA_real_)
        }
      }
    }
  )

  if (all(sapply(parm_gas, is.na))) {
    starts <- parm_gas
    par_DF$start <- starts[par_DF$param_name]
    return(par_DF)
  }

  # Set parameters needed for model
  Fup <- parm_gas[["Funbound.plasma"]]
  Rblood2plasma <- parm_gas[["Rblood2plasma"]]
  Clint <- parm_gas[["Clint"]]
  Kblood2air <- parm_gas[["Kblood2air"]]
  Q_totli <- (parm_gas[["Qliverf"]] + parm_gas[["Qgutf"]]) * parm_gas[["Qcardiacc"]]
  Q_gfr <- parm_gas[["Qgfrc"]]
  Q_alv <- parm_gas[["Qalvc"]]
  kgutabs <- parm_gas[["kgutabs"]]
  Fgutabs <- parm_gas[["Fabsgut"]]
  BW <- parm_gas[["BW"]]

  liver_data <- httk::tissue.data[
    httk::tissue.data$Tissue == "liver" &
      tolower(httk::tissue.data$Species) == this_species,
  ]

  liver_mass <- liver_data[liver_data$variable == "Density (g/cm^3)", "value"] *
    liver_data[liver_data$variable == "Vol (L/kg)", "value"] * 1E3 # mL to L

  # Get starting Concs from data

  # Work only with detects for these rough estimates
  tmpdat <- subset(data,
                   Detect %in% TRUE)

  # Split into IV and PO
  ivdat <- subset(tmpdat,
                  Route %in% "iv")
  podat <- subset(tmpdat,
                  Route %in% "oral")

  # Set a Fup specific to the liver for clearance
  if (!restrictive) {
    Fup_hep <- 1
  } else {
    Fup_hep <- Fup
  }

  # Convert L/h/kg BW^(3/4) to L/h
  Q_totli_2 <- Q_totli / (BW^(3/4))
  Q_gfr_2 <- Q_gfr / (BW^(3/4))
  Q_alv_2 <- Q_alv / (BW^(3/4))

  Clint_hep <- Clint * (6.6 * BW * liver_mass) / 1E6 # Convert to L/hr
  Clhep <- Q_totli_2 * (Fup_hep * Clint_hep / Rblood2plasma) /
    (Q_totli_2 + (Fup_hep * Clint_hep / Rblood2plasma))
  # Need to include Fup for renal clearance
  Clren <- Fup * Q_gfr_2 / Rblood2plasma
  Clair <- (Rblood2plasma * Q_alv_2 / Kblood2air)

  Cltot <- Clren + Clhep + Clair


  # Quick and dirty:
  # IV data estimates, if IV data exist
  if (nrow(ivdat) > 0) {

    # assume that midpoint of time is one half-life, so kelim = log(2)/(midpoint of time).
    halflife <- mean(range(ivdat$Time))
    kelim <- log(2) / halflife

    # Vdist: calculate this based on Cltot/kelim
    Vdist <- Cltot / kelim
  }

  if (nrow(podat) > 0) {
    # if PO data exist, then we can get ka and Fgutabs/Vdist
    # get peak time
    tCmax <- get_peak(x = podat$Time,
                      y = log10(podat$Conc / podat$Dose))
    tmax <- tCmax[[1]]
    Cmax <- tCmax[[2]]


    # if no IV data, then calculate kelim from oral data
    if (nrow(ivdat) == 0) {
      # and assume that midpoint of time is one half-life, so kelim = log(2)/(midpoint of time).
      halflife <- mean(range(podat$Time))
      kelim <- log(2) / halflife
      Vdist <- Cltot / kelim
    }
  }


  starts <- c("Q_totli" = Q_totli,
              "Q_gfr" = Q_gfr,
              "Q_alv" = Q_alv,
              "Kblood2air" = Kblood2air,
              "BW" = BW,
              "Fup" = Fup,
              "Clint" = Clint,
              "kgutabs" = kgutabs,
              "Vdist" = Vdist,
              "Fgutabs" = Fgutabs,
              "liver_mass" = liver_mass,
              "Rblood2plasma" = Rblood2plasma)

  par_DF$start <- starts[par_DF$param_name]

  return(par_DF)
}
