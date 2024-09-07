#' Get starting values for flat model
#'
#' Derive starting values for flat model parameters from available data
#'
#' This function is called internally by [get_params_1comp()] and should
#' generally not be called directly by the user.
#'
#'The full set of model parameters for the flat model includes
#'`Vdist`,`Fgutabs`, and `Rblood2plasma`. Whether each one can be estimated from
#'the data depends on what routes of administration are included in the data.
#'
#'The numerical optimizer requires starting guesses for the value of each
#'parameter to be estimated from the data. Default starting guesses are derived from the available data.
#'
#'These are intended to be *very* rough starting guesses, so the algorithm here
#'is extremely naive. This function is not itself intended to produce valid
#'estimates for any of the model parameters, and it is highly unlikely to do so.
#'
#'The derivation process is as follows.
#'
#'First, data are filtered to exclude any non-detects.
#'
#'Then, data are split by route of administration, into an IV data set and an oral data
#'set. (It is possible that either IV or oral data may not be
#'available for a chemical.)
#'
#'If IV data exist, then only IV data are used to derive a starting estimate for
#'`Vdist`. Concentrations are dose-normalized (divided by their corresponding
#'dose) and log10-transformed. The mean dose-normalized, log10-transformed
#'concentration is calculated (call it `Cmean_log10`). `Vdist` starting value is
#'then derived as `1/(10^Cmean_log10)` .
#'
#'If any oral data exist (whether or not IV data also exist), then the oral data
#'are used to derive a starting value for `Fgutabs_Vdist`. Concentrations are dose-normalized (divided by their corresponding
#'dose) and log10-transformed. The mean dose-normalized, log10-transformed
#'concentration is calculated (call it `Cmean_log10`). `Fgutabs_Vdist` starting value is
#'then set equal to `10^Cmean_log10` .
#'
#'#Starting value for `Rblood2plasma`
#'
#' If both blood and plasma data are available, then the starting value for `Rblood2plasma` is derived as follows.
#'
#' If IV data are available for both blood and plasma, then the starting value
#' for `Rblood2plasma` is derived as the ratio of `Vdist` for blood data and
#' `Vdist` for plasma data.
#'
#' If oral data, but not IV data, are available for both blood and plasma, then
#' the starting value for `Rblood2plasma` is derived as the ratio of
#' `Fgutabs_Vdist` for plasma data and `Fgutabs_Vdist` for blood data.
#'
#' If only blood data or only plasma data are available, then the starting value for `Rblood2plasma` is set at a constant 1.
#'
#'@param data The data set to be fitted (e.g. the result of [preprocess_data()])
#' @param par_DF A `data.frame` with the following variables (e.g., as produced by [get_params_flat()])
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
#' @family flat model functions
#' @family get_starts functions
#' @family built-in model functions
#'
get_starts_flat <- function(data,
                            par_DF){

  Vdist <- NA_real_
  Fgutabs_Vdist <- NA_real_
  Fgutabs <- NA_real_
  Rblood2plasma <- 1

  #Work only with detects for these rough estimates
  tmpdat <- subset(data,
                   Detect %in% TRUE)

  #Split into IV and PO
  ivdat <- subset(tmpdat,
                  Route %in% "iv")
  podat <- subset(tmpdat,
                  Route %in% "oral")

  has_iv <- any(tmpdat$Route %in% "iv")
  has_po <- any(tmpdat$Route %in% "oral")
  has_plasma <- any(tmpdat$Media %in% "plasma")
  has_blood <- any(tmpdat$Media %in% "blood")
  has_iv_plasma <- any(tmpdat$Media %in% "plasma" &
                      tmpdat$Route %in% "iv")
  has_iv_blood <- any(tmpdat$Media %in% "blood" &
                         tmpdat$Route %in% "iv")
  has_po_plasma <- any(tmpdat$Media %in% "plasma" &
                         tmpdat$Route %in% "oral")
  has_po_blood <- any(tmpdat$Media %in% "blood" &
                        tmpdat$Route %in% "oral")

  Vdist_plasma_log10 <- NA_real_
  Vdist_blood_log10 <- NA_real_

  Fgutabs_Vdist_plasma_log10 <- NA_real_
  Fgutabs_Vdist_blood_log10 <- NA_real_

  Rblood2plasma_iv_log10 <- NA_real_
  Rblood2plasma_po_log10 <- NA_real_

  if (has_iv %in% TRUE) {
    if (has_iv_plasma %in% TRUE) {
      #get Vdist for plasma
      Cmean_plasma_log10 <- with(subset(ivdat, Media %in% "plasma"),
                                 mean(log10(Conc/Dose), na.rm = TRUE))
      Vdist_plasma_log10 <- -Cmean_plasma_log10
    }
    if (has_iv_blood %in% TRUE) { #if both blood and plasma IV data, estimate Rblood2plasma
      #get Vdist for blood

      Cmean_blood_log10 <- with(subset(ivdat,
                                       Media %in% "blood"),
                                mean(log10(Conc/Dose), na.rm = TRUE))
      Vdist_blood_log10 <- -Cmean_blood_log10
    }
    if (has_iv_plasma %in% TRUE & has_iv_blood %in% TRUE) {
      Rblood2plasma_iv_log10 <- Vdist_plasma_log10 - Vdist_blood_log10
    }
  }

  if (has_po %in% TRUE) {
    if (has_po_plasma %in% TRUE) {
      #get Fgutabs_Vdist for plasma
      Fgutabs_Vdist_plasma_log10 <- with(subset(podat,Media %in% "plasma"),
                                         mean(log10(Conc/Dose),na.rm = TRUE))
    }
    if (has_po_blood %in% TRUE) {
      #get Fgutabs_Vdist for blood
      Fgutabs_Vdist_blood_log10 <- with(subset(podat, Media %in% "blood"),
                                        mean(log10(Conc/Dose), na.rm = TRUE))
    }
    if (has_po_plasma %in% TRUE & has_po_blood %in% TRUE) {
      Rblood2plasma_po_log10 <- Fgutabs_Vdist_blood_log10 + Fgutabs_Vdist_blood_log10
    }
  }

  if (has_iv %in% TRUE) {
    Vdist <- 10^(mean(c(Vdist_plasma_log10, Vdist_blood_log10),
                      na.rm = TRUE))
  }

  if (has_po %in% TRUE) {
    Fgutabs_Vdist <- 10^(mean(c(Fgutabs_Vdist_plasma_log10, Fgutabs_Vdist_blood_log10),
                              na.rm = TRUE))
  }

  if (has_iv %in% TRUE & has_po %in% TRUE) {
    Fgutabs <- Fgutabs_Vdist * Vdist
  }

  if (has_plasma %in% TRUE & has_blood %in% TRUE) {
    Rblood2plasma <- 10^(mean(c(Rblood2plasma_iv_log10, Rblood2plasma_po_log10),
                              na.rm = TRUE))
  }

  #update starting Concs
  par_DF["Vdist", "start"] <- Vdist
  par_DF["Fgutabs_Vdist", "start"] <- Fgutabs_Vdist
  par_DF["Fgutabs", "start"] <- Fgutabs
  par_DF["Rblood2plasma", "start"] <- Rblood2plasma

  return(par_DF)
}
