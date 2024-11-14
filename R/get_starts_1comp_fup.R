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
#' # Starting value for `kelim`
#'
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
#' # Starting value for `Vdist`
#'
#' Using a calculated value for total clearance, `Cl_tot`, `Vdist` is estimated
#' by dividing this by the estimation of `kelim`.
#'
#' # Starting value for `kgutabs`
#'
#' If oral data exist (whether or not IV data also exist), then the oral data
#' are used to derive a starting value for `kgutabs`.
#'
#' First, concentrations are dose-normalized by dividing them by their corresponding
#' doses. Then the normalized concentrations are log10-transformed.
#'
#' The time of peak concentration (`tmax`), and the median (normalized,
#' log-transformed) peak concentration (`Cmax_log10`), are identified using [get_peak()].
#'
#' As a very rough guess,`tmax` is assumed to
#' occur at one absorption half-life. Under this assumption, `kgutabs` is equal
#' to `log(2)/tmax`, and this is taken as the starting value.
#'
#' # Starting value for `Fgutabs_Vdist`
#'
#' If any oral data exist (whether or not IV data also exist), then the oral data
#' are used to derive a starting value for `Fgutabs_Vdist`.
#'
#' If the kinetics obey a one-compartment model, then if concentrations are
#' dose-normalized, log-transformed, and plotted vs. time, then at late time
#' points (after concentration has peaked), the concentration vs. time
#' relationship will approach a straight line with slope `-kelim`.
#'
#' If this straight line is extrapolated back to time 0, then the resulting
#' intercept (call it `A`), expressed on the natural scale, is equal to
#' `Fgutabs_Vdist * kgutabs/(kgutabs-kelim)`. See
#' https://www.boomer.org/c/p4/c09/c0902.php .
#'
#' Roughly, we approximate `A` on the log10 scale by extrapolating back from the peak along a
#' straight line with slope `-kelim`, using the previously-derived starting
#' value for `kelim`. So `log10(A) = Cmax_log10 + kelim*tmax`.
#'
#' Using the previously-derived starting values for `kgutabs` and `kelim`, then,
#' the starting value for `Fgutabs_Vdist` can be derived as `A * (kgutabs-kelim)/kgutabs`.
#'
#' # Starting value for `Fgutabs`
#'
#' If both oral and IV data exist, then the derived starting values for `Vdist`
#' (from the IV data) and `Fgutabs_Vdist` (from the oral data) are multiplied to
#' yield a derived starting value for `Fgutabs`.
#'
#' #Starting value for `Rblood2plasma`
#'
#' The starting value for `Rblood2plasma` is set to the value given by [httk::parameterize_gas_pbtk()].
#'
#' @param data The data set to be fitted (e.g. the result of [preprocess_data()])
#' @param par_DF A `data.frame` with the following variables (e.g., as produced by [get_params_1comp()])
#' - `param_name`: Character: Names of the model parameters
#' - `param_units`: Character: Units of the model parameters
#' - `optimize_param`: TRUE if each parameter is to be estimated from the data; FALSE otherwise
#' - `use_param`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
#' -`lower_bounds`: Numeric: The lower bounds for each parameter
#' - `upper_bounds`: Numeric: The upper bounds for each parameter
#'
#' @return The same `data.frame` as `par_DF`, with an additional variable
#'  `starts` containing the derived starting value for each parameter. If a
#'  parameter cannot be estimated from the available data, then its starting value
#'  will be `NA_real_`
#' @import httk
#' @author Caroline Ring
#' @family 1-compartment model functions
#' @family get_starts functions
#' @family built-in model functions

get_starts_1comp_fup <- function(data,
                                 par_DF) {
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

  parm_gas <- suppressMessages(
    suppressWarnings(
      httk::parameterize_gas_pbtk(
        dtxsid = unique(data[["Chemical"]]),
        species = this_species,
        restrictive.clearance = TRUE)))

  Fup <- parm_gas[["Funbound.plasma"]]

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


  starts <- c("Q_totli" = Q_totli,
              "Q_gfr" = Q_gfr,
              "Fup" = Fup,
              "Clint" = Clint,
              "kgutabs" = kgutabs,
              "Vdist" = Vdist,
              "Fgutabs_Vdist" = Fgutabs_Vdist,
              "Fgutabs" = Fgutabs,
              "Rblood2plasma" = Rblood2plasma)

  par_DF$start <- starts[par_DF$param_name]

  return(par_DF)
}
