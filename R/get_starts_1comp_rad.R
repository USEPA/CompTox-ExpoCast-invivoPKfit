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
  Fgutabs_Vdist <- NA_real_
  Fup <- 1

  Q_gfr <- httk::physiology.data %>%
    dplyr::filter(Parameter %in% "GFR") %>%
    tidyr::pivot_longer(cols = Mouse:Monkey,
                        names_to = "Species",
                        values_to = "param_value")
  Q_gfr <- setNames(object = Q_gfr[["param_value"]],
                    nm = tolower(Q_gfr[["Species"]]))

  Q_totli <- httk::tissue.data %>%
    dplyr::filter(variable %in% "Flow (mL/min/kg^(3/4))",
                  Tissue %in% "liver")
  Q_totli <- setNames(object = Q_totli[["value"]],
                      nm = tolower(Q_totli[["Species"]]))

  names_Q_gfr <- names(Q_gfr)
  this_species <- unique(data$Species)

  if (this_species %in% names_Q_gfr) {
    Q_gfr <- Q_gfr[[this_species]] * (60 / 1000) # Assumes L/h/kg are standard units
    Q_totli <- Q_totli[[this_species]] * (60 / 1000)
  } else {
    Q_gfr <- Q_gfr[["human"]] * (60 / 1000)
    Q_totli <- Q_totli[["human"]] * (60 / 1000)
    message("Species not in database, using human values for Q_gfr & Q_totli")

  }

  parm_gas <- tryCatch(expr = {
    suppressMessages(
      suppressWarnings(
        httk::parameterize_3comp2(
          dtxsid = unique(data[["Chemical"]]),
          species = this_species,
          default.to.human = TRUE,
          restrictive.clearance = restrictive)))
  }, error = function(e) {
    warning("There are no parameters for: ", unique(data[["Chemical"]]))
    response = tolower(
      trimws(
        readline(prompt = "Do you wish to substitute Bis-phenol A's parameters?")
      )
    )
    if (response %in% c("y", "yes")) {
      suppressMessages(
        suppressWarnings(
          httk::parameterize_3comp2(
            "bisphenol a",
            species = this_species,
            default.to.human = TRUE,
            restrictive.clearance = restrictive)))
    } else {
      stop("You have elected not to substitute, prefit cannot continue.")
    }
  }
  )


  if (restrictive) {
    Fup <- parm_gas[["Funbound.plasma"]]
  }

  Rblood2plasma <- parm_gas[["Rblood2plasma"]]

  Clint <- parm_gas[["Clint"]]

  Kblood2air <- parm_gas[["Kblood2air"]]

  Qalvc <- parm_gas[["Qalvc"]]


  # Get starting Concs from data

  # Work only with detects for these rough estimates
  tmpdat <- subset(data,
                   Detect %in% TRUE)

  # Split into IV and PO
  ivdat <- subset(tmpdat,
                  Route %in% "iv")
  podat <- subset(tmpdat,
                  Route %in% "oral")

  Cl_hep <- Q_totli * Fup * Clint / (Q_totli + (Fup * Clint / Rblood2plasma))
  Cl_tot <- Q_gfr + Cl_hep + (Rblood2plasma * Qalvc / Kblood2air)


  # Quick and dirty:
  # IV data estimates, if IV data exist
  if (nrow(ivdat) > 0) {

    # assume that midpoint of time is one half-life, so kelim = log(2)/(midpoint of time).
    halflife <- mean(range(ivdat$Time))
    kelim <- log(2) / halflife

    # Vdist: calculate this based on Cltot/kelim
    Vdist <- Cl_tot / kelim
  }

  if (nrow(podat) > 0) {
    # if PO data exist, then we can get ka and Fgutabs/Vdist
    # get peak time
    tCmax <- get_peak(x = podat$Time,
                      y = log10(podat$Conc / podat$Dose))
    tmax <- tCmax[[1]]
    Cmax <- tCmax[[2]]

    # assume peak time occurs at 1 absorption halflife
    # so kgutabs = log(2)/tmax
    kgutabs <- log(2) / tmax

    # if no IV data, then calculate kelim from oral data
    if (nrow(ivdat) == 0) {
      # and assume that midpoint of time is one half-life, so kelim = log(2)/(midpoint of time).
      halflife <- mean(range(podat$Time))
      kelim <- log(2) / halflife
      Vdist <- Cl_tot / kelim
    }


    # then extrapolate back from Cmax to time 0 with slope -kelim
    Fgutabs_Vdist <- 10^((Cmax + kelim * tmax)) * (kgutabs - kelim) / (kgutabs)
    # if we had IV data, then we had a Vdist estimate, so we can estimate Fgutabs too
    Fgutabs <- Fgutabs_Vdist * Vdist
  }

  # Get the initial guess for fraction recovered
  # Does the excreta data have any non-detects?
  # Get the time at which the maximum cumulative concentration occurs
  data_exc <- subset(data, subset = Media %in% 'excreta')
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
