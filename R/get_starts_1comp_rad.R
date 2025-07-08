#' Get starting values for 1-compartment model with specific clearance
#'
#' Derive starting values for 1-compartment model parameters from available data
#'
#' This function is called internally by [get_params_1comp()] and should
#' generally not be called directly by the user.
#'
#' The full set of model parameters for the 1-compartment model includes `Vdist`,
#' `kelim`, `kgutabs`, `Fgutabs`, and `Rblood2plasma`. Whether each one can be
#' estimated from the data depends on what routes of administration are included
#' in the data. However, in this version of the one-compartment model,
#' we use the liver and glomerular flow rates (`Q_totli` and `Q_gfr`, respectively)
#' provided by [httk] as well as intrisic clearance `Clint` and fraction unbound
#' in plasma, `Fup`. These starting values are established and then are used
#' to calculate the other parameters in the 1-compartment model.
#'
#' The numerical optimizer requires starting guesses for the value of each
#' parameter to be estimated from the data. Default starting guesses are derived from the available data.
#'
#' These are intended to be *very* rough starting guesses, so the algorithm here
#' is extremely naive. This function is not itself intended to produce valid
#' estimates for any of the model parameters, and it is highly unlikely to do so.
#'
#' The derivation process is as follows.
#'
#' First, data are filtered to exclude any non-detects.
#'
#' Then, data are split by route of administration, into an IV data set and an oral data
#' set. (It is possible that either IV or oral data may not be
#' available for a chemical.)
#'
#' @inheritSection get_starts_1comp_cl Starting value for `Vdist`
#' @section Starting value for `Frec`:
#'
#' The empirical value of `Frec` is attainable when the experimental collection
#' of excreta goes on until there is no detectable compound. In this case,
#' the values will be equal to the maximum cumulative amount recovered divided
#' by `Dose` or `Dose` * `Fgutabs` for intravenous and oral administration, respectively.
#'
#' If all `excreta` data are 'detectable', the assumption that we can necessarily
#' obtain an empirical estimate of `Frec` may not hold. While in this case still,
#' the maximum is used, the model fitting steps will help realize a better estimate
#' for this fraction recovered.
#'
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
#' @author Gilberto Padilla Mercado, Caroline Ring
#' @family 1-compartment radiation model functions
#' @family get_starts functions
#' @family built-in model functions
#'
get_starts_1comp_rad <- function(data,
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
  Fup <- 1
  Frec <- NA_real_
  Fgutabs_Vdist <- NA_real_

  Q_gfr <- httk::physiology.data[
    httk::physiology.data$Parameter == "GFR",
    c("Mouse", "Rat", "Dog", "Human", "Rabbit", "Monkey")
  ]
  Q_gfr <- setNames(object = unlist(Q_gfr),
                    nm = tolower(names(Q_gfr)))


  Q_totli <- httk::tissue.data[
    httk::tissue.data$Tissue == "liver" &
      httk::tissue.data$variable == "Flow (mL/min/kg^(3/4))",
    c("Species", "value")
  ]
  Q_totli <- setNames(object = Q_totli[["value"]],
                      nm = tolower(Q_totli[["Species"]]))

  Q_alv <- httk::physiology.data[
    httk::physiology.data$Parameter == "Pulmonary Ventilation Rate",
    c("Mouse", "Rat", "Dog", "Human", "Rabbit", "Monkey")
  ]
  Q_alv <- setNames(object = unlist(Q_alv),
                    nm = tolower(names(Q_alv)))

  # Get species-specific flow rates, or default to human
  names_Q_gfr <- names(Q_gfr)
  names_Q_alv <- names(Q_alv)
  this_species <- unique(data$Species)

  if (this_species %in% names_Q_gfr) {
    Q_gfr <- Q_gfr[[this_species]] * (60 / 1000) # Assumes L/h/kg are standard units
    Q_totli <- Q_totli[[this_species]] * (60 / 1000)
  } else {
    Q_gfr <- Q_gfr[["human"]] * (60 / 1000)
    Q_totli <- Q_totli[["human"]] * (60 / 1000)
    message("Species not in database, using human values for Q_gfr & Q_totli")
  }

  if (this_species %in% names_Q_alv) {
    Q_alv <- Q_alv[[this_species]]
  } else {
    Q_alv <- Q_alv[["human"]]
    message("Species not in database, using human values for Q_alv")
  }

  parm_gas <- tryCatch(
    expr = {
      httk::parameterize_3comp2(
        dtxsid = unique(data[["Chemical"]]),
        species = this_species,
        default.to.human = TRUE,
        restrictive.clearance = restrictive
      ) |>
        suppressWarnings() |>
        suppressMessages()
    }, error = function(e) {
      c("Q_totli" = NA_real_,
        "Q_gfr" = NA_real_,
        "Q_alv" = NA_real_,
        "Kblood2air" = NA_real_,
        "Fup" = NA_real_,
        "Clint" = NA_real_,
        "kgutabs" = NA_real_,
        "Vdist" = NA_real_,
        "Fgutabs_Vdist" = NA_real_,
        "Fgutabs" = NA_real_,
        "Rblood2plasma" = NA_real_,
        "Frec" = NA_real_)
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
  Q_alv <- parm_gas[["Qalvc"]]
  kgutabs <- parm_gas[["kgutabs"]]
  Fgutabs <- parm_gas[["Fabsgut"]]


  # Get starting Concs from data

  # Work only with detects for these rough estimates
  tmpdat <- subset(data,
                   Detect %in% TRUE)

  # Split into IV and PO
  ivdat <- subset(tmpdat,
                  Route %in% "iv" & Media %in% c("blood", "plasma"))
  podat <- subset(tmpdat,
                  Route %in% "oral" & Media %in% c("blood", "plasma"))

  # Set a Fup specific to the liver for clearance
  if (!restrictive) {
    Fup_hep <- 1
  } else {
    Fup_hep <- Fup
  }
  Clhep <- Q_totli * Fup_hep * Clint / (Q_totli + (Fup_hep * Clint / Rblood2plasma))
  # Need to include Fup for renal clearance
  Clren <- Fup * Q_gfr
  Clair <- (Rblood2plasma * Q_alv / Kblood2air)

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

  # Get the initial guess for fraction recovered
  # Does the excreta data have any non-detects?
  # Get the time at which the maximum cumulative concentration occurs
  data_exc <- subset(data, subset = Media %in% "excreta")
  tmax_exc <- data_exc[
    which(data_exc[["Conc"]] == max(data_exc[["Conc"]])),
    "Time"
    ]

  # Take the median of the concentrations at that timepoint
  # (accounting for multiple individual observations)
  Frec <- median(data_exc[data_exc[["Time"]] == tmax_exc, "Conc"],
                 na.rm = TRUE)

  starts <- c("Q_totli" = Q_totli,
              "Q_gfr" = Q_gfr,
              "Fup" = Fup,
              "Clint" = Clint,
              "kgutabs" = kgutabs,
              "Vdist" = Vdist,
              "Fgutabs_Vdist" = Fgutabs_Vdist,
              "Fgutabs" = Fgutabs,
              "Rblood2plasma" = Rblood2plasma,
              "Frec" = Frec)

  par_DF$start <- starts[par_DF$param_name]

  return(par_DF)
}
