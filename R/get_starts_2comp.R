#' Get starting values for 2-compartment model
#'
#' Derive starting values for 2-compartment model parameters from available data
#'
#' This function is called internally by [get_params_2comp()] and should
#' generally not be called directly by the user.
#'
#' The full set of model parameters for the 2-compartment model includes `V1`,
#' `kelim`, `k12`, `k21`, `kgutabs`, and `Fgutabs`. Whether each one can be
#'estimated from the data depends on what routes of administration are included
#'in the data.
#'
#'The numerical optimizer requires starting guesses for the value of each
#'parameter to be estimated from the data. Default starting guesses are derived from the available data.
#'
#'These are intended to be *very* rough starting guesses, so the algorithm here
#'is extremely naive. This function is not itself intended to produce valid
#'estimates for any of the model parameters, and it is highly unlikely to do so.
#'
#'At present, the starting guesses for the 2-compartment model are derived in
#'the same way as for the 1-compartment model, for the parameters that are
#'common to both. That is, the data are assumed to obey a 1-compartment model to
#'derive starting guesses for `kelim`, `V1`, `kgutabs`, `Fgutabs_V1`, and
#'`Fgutabs`.
#'
#'Then, starting values for `k12` and `k21` are arbitrarily set to
#'0.1 and 0.5, respectively.
#'
#'The following description of the derivation process is therefore identical to
#'that for [get_starts_1comp()].
#'
#'The derivation process
#'
#'First, data are filtered to exclude any non-detects.
#'
#'Then, data are split by route of administration, into an IV data set and an oral data
#'set. (It is possible that either IV or oral data may not be
#'available for a chemical.)
#'
#'# Starting value for `kelim`
#'
#'If IV data exist, then only IV data are used to derive starting estimates for
#'`kelim`, even if oral data also exist.
#'
#'If only oral data exist, then the oral data are used to derive a starting
#'estimate for `kelim`.
#'
#'Whichever data set is used (IV or oral), the starting value for `kelim` is
#'derived by assuming that the range of observed time values in the data set
#'spans two elimination half-lives. This implies that the elimination half-life
#'is equal to the midpoint of observed time values, and that the starting value
#'for the elimination time constant `kelim` is therefore `log(2)` divided by the
#'midpoint of observed time values.
#'
#'Of course, this assumption is unlikely to be correct. However, we hope that it
#'will yield a starting guess for `kelim` that is at least on the right order of
#'magnitude.
#'
#' # Starting value for `V1`
#'
#' If IV data exist, then only IV data are used to derive a starting estimate
#' for `V1`.
#'
#' This starting estimate is derived by assuming that the IV data obey a
#' one-compartment model, which means that when concentrations are
#' dose-normalized and log10-transformed and plotted against time, they will
#' follow a straight line with slope `-kelim`.
#'
#' First, concentrations are dose-normalized by dividing them by their corresponding
#' doses. Then the normalized concentrations are log10-transformed.
#'
#' From all observations at the earliest observed time point in the data set
#' (call it `tmin`), the median of the dose-normalized, log10-transformed
#' concentrations is calculated; call it `C_tmin`. (The median is used, rather
#' than the mean, in an attempt to be more robust to outliers.)
#'
#' If the earliest observed time point is not at time = 0, then the
#' dose-normalized, log10-transformed concentration at time = 0 is extrapolated
#' by drawing a straight line with slope `-kelim` back from `C_tmin`, where the
#' value of `kelim` is the starting value derived as in the previous section.
#'
#' This extrapolated concentration at time t = 0 is called `A_log10`. `A_log10`
#' represents the expected body concentration immediately after IV injection of
#' a unit dose (under the assumption that TK obeys a one-compartment model).
#'
#' Then, the volume of distribution `V1` is derived as `1/(10^A_log10)`. In
#' other words, `V1` is the volume that would be required to produce a
#' concentration equal to `A_log10` after injecting a unit dose.
#'
#' (No starting value for `V1` can be derived with only oral data,
#'but none is needed, because with only oral data, `V1` will not be estimated
#'from the data).
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
#' # Starting value for `Fgutabs_V1`
#'
#'If any oral data exist (whether or not IV data also exist), then the oral data
#'are used to derive a starting value for `Fgutabs_V1`.
#'
#' If the kinetics obey a one-compartment model, then if concentrations are
#' dose-normalized, log-transformed, and plotted vs. time, then at late time
#' points (after concentration has peaked), the concentration vs. time
#' relationship will approach a straight line with slope `-kelim`.
#'
#' If this straight line is extrapolated back to time 0, then the resulting
#' intercept (call it `A`), expressed on the natural scale, is equal to
#' `Fgutabs_V1 * kgutabs/(kgutabs-kelim)`. See
#' https://www.boomer.org/c/p4/c09/c0902.php .
#'
#' Roughly, we approximate `A` on the log10 scale by extrapolating back from the peak along a
#' straight line with slope `-kelim`, using the previously-derived starting
#' value for `kelim`. So `log10(A) = Cmax_log10 + kelim*tmax`.
#'
#' Using the previously-derived starting values for `kgutabs` and `kelim`, then,
#' the starting value for `Fgutabs_V1` can be derived as `A * (kgutabs-kelim)/kgutabs`.
#'
#'# Starting value for `Fgutabs`
#'
#'If both oral and IV data exist, then the derived starting values for `V1`
#'(from the IV data) and `Fgutabs_V1` (from the oral data) are multiplied to
#'yield a derived starting value for `Fgutabs`.
#'#Starting value for `Rblood2plasma`
#'
#'The starting value for `Rblood2plasma` is always set at a constant 1.
#'
#'@param data The data set to be fitted (e.g. the result of [preprocess_data()])
#' @param par_DF A `data.frame` with the following variables (e.g., as produced by [get_params_2comp()])
#' - `param_name`: Character: Names of the model parameters
#' - `param_units`: Character: Units of the model parameters
#' - `optimize_param`: TRUE if each parameter is to be estimated from the data; FALSE otherwise
#' - `use_param`: TRUE if each parameter is to be used in evaluating the model; FALSE otherwise
#' -`lower_bounds`: Numeric: The lower bounds for each parameter
#' - `upper_bounds`: Numeric: The upper bounds for each parameter
#'
#'@return The same `data.frame` as `par_DF`, with an additional variable
#'  `starts` containing the derived starting value for each parameter. If a
#'  parameter cannot be estimated from the available data, then its starting value
#'  will be `NA_real_`
#'
#' @author Caroline Ring
#' @family 2-compartment model functions
#' @family get_starts functions
#' @family built-in model functions
get_starts_2comp <- function(data,
                             par_DF){

  kelim <- NA_real_
  kgutabs <- NA_real_
  V1 <- NA_real_
  Fgutabs_V1 <- NA_real_
  k12 <- NA_real_
  k21 <- NA_real_
  Fgutabs <- NA_real_
  Rblood2plasma <- 1

  # Get starting Concs from data

  # Get starting Concs from data

  #Work only with detects for these rough estimates
  tmpdat <- subset(data,
                   Detect %in% TRUE)

  #Split into IV and PO
  ivdat <- subset(tmpdat,
                  Route %in% "iv")
  podat <- subset(tmpdat,
                  Route %in% "oral")

  # Quick and dirty:
  #IV data estimates, if IV data exist
  if(nrow(ivdat)>0){

    #assume that midpoint of time is one half-life, so kelim = log(2)/(midpoint of time).
    halflife <- mean(range(ivdat$Time))
    kelim <- log(2)/halflife

    #V1: extrapolate back from conc at min time at a slope of -kelim to get the intercept
    #then V1 = 1/intercept
    C_tmin <- with(subset(ivdat, Time == min(Time)),
                   median(log10(Conc/Dose)))
    A_log10 <- C_tmin + kelim*min(ivdat$Time)
    V1 <- 1/(10^A_log10)
  }

  if(nrow(podat)>0){
    #if PO data exist, then we can get ka and Fgutabs/V1
    #get peak time
    tCmax <- get_peak(x = podat$Time,
                      y = log10(podat$Conc/podat$Dose))
    tmax <- tCmax[[1]]
    Cmax <- tCmax[[2]]

    #assume peak time occurs at 1 absorption halflife
    #so kgutabs = log(2)/tmax
    kgutabs <- log(2)/tmax

    #if no IV data, then calculate kelim from oral data
    if(nrow(ivdat)==0){
      #and assume that midpoint of time is one half-life, so kelim = log(2)/(midpoint of time).
      halflife <- mean(range(podat$Time))
      kelim <- log(2)/halflife
    }

    #then extrapolate back from Cmax to time 0 with slope -kelim
    Fgutabs_V1 <- 10^((Cmax + kelim*tmax))*(kgutabs - kelim)/(kgutabs)

    if(nrow(ivdat)>0){
      #if we had IV data, then we had a V1 estimate, so we can estimate Fgutabs too
      Fgutabs <- Fgutabs_V1 * V1
    }
  }

  #arbitrary until I implement something more sophisticated
  k21 <- 0.1
  k12 <- 0.5

  starts <- c("kelim" = kelim,
              "kgutabs" = kgutabs,
              "k12" = k12,
              "k21" = k21,
              "V1" = V1,
              "Fgutabs_V1" = Fgutabs_V1,
              "Fgutabs" = Fgutabs,
              "Rblood2plasma" = Rblood2plasma)

par_DF$start <- starts[par_DF$param_name]

  return(par_DF)
}
