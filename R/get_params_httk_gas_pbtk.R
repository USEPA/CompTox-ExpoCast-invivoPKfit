#' Get parameters to fit `httk`'s `gas_pbtk` PBPK model
#'
#' Get parameters to fit the `gas_pbtk`
#' model from the `httk` package (Wambaugh, Schacht, and Ring. 2025).
#'
#' @section Required parameters:
#' These are given by the \link[httk]{parameterize_3comp2} function in `httk`.
#' Furthermore, they are transformed to a vector during hte prefitting process.
#'
#' @param data The data set to be fitted (e.g. the result of [preprocess_data()])
#' @param lower_bound A mapping specified using a call to [alist()],
#'  giving the lower bounds for each variable, as expressions which may include
#'  variables in `data`.
#' @param upper_bound A mapping specified using a call to [alist()],
#'  giving the upper bounds for each variable, as expressions which may include
#'  variables in `data`.
#' @param param_units A mapping specified using a call to [alist()],
#'  giving the units for each variable, as expressions which may include
#'  variables in `data`.
#' @param restrictive A logical value (TRUE or FALSE. Default: FALSE) that says whether the
#' assumption is that the clearance is restrictive or non-restrictive
#'
#' @return A vector of blood or plasma concentration values  corresponding
#'  to `time`.
#' @author Gilberto Padilla Mercado
#' @export get_params_httk_gas_pbtk
#' @family built-in model functions
#' @family httk model functions
#' @family model concentration functions
#'
get_params_httk_gas_pbtk <- function(
    data,
    lower_bound = alist(
      BW = NA,
      Caco2.Pab = NA,
      Caco2.Pab.dist = NA,
      Clint = 0,
      Clint.dist = NA,
      Clmetabolismc = NA,
      Funbound.plasma = 0,
      Funbound.plasma.dist = NA,
      Funbound.plasma.adjustment = NA,
      Fabsgut = NA,
      Fhep.assay.correction = NA,
      hematocrit = NA,
      Kgut2pu = NA,
      Krbc2pu = NA,
      kgutabs = NA,
      Kkidney2pu = NA,
      Klung2pu = NA,
      km = NA,
      Kmuc2air = NA,
      Kliver2pu = NA,
      Krest2pu = NA,
      Kblood2air = NA,
      kUrtc = NA,
      liver.density = NA,
      logHenry = NA,
      million.cells.per.gliver = NA,
      MW = NA,
      Pow = NA,
      pKa_Donor = NA,
      pKa_Accept = NA,
      MA = NA,
      Qcardiacc = NA,
      Qgfrc = NA,
      Qgutf = NA,
      Qliverf = NA,
      Qalvc = NA,
      Qkidneyf = NA,
      Qlungf = NA,
      Rblood2plasma = NA,
      Vgutc = NA,
      Vliverc = NA,
      Vartc = NA,
      Vkidneyc = NA,
      Vlungc = NA,
      vmax = NA,
      Vmucc = NA,
      Vvenc = NA,
      Vrestc = NA,
      KFsummary = NA,
      Fprotein.plasma = NA,
      fabs.oral = NA,
      Qgut_ = NA,
      Qintesttransport = NA
    ),
    upper_bound = alist(
      BW = NA,
      Caco2.Pab = NA,
      Caco2.Pab.dist = NA,
      Clint = 1000,
      Clint.dist = NA,
      Clmetabolismc = NA,
      Funbound.plasma = 1,
      Funbound.plasma.dist = NA,
      Funbound.plasma.adjustment = NA,
      Fabsgut = NA,
      Fhep.assay.correction = NA,
      hematocrit = NA,
      Kgut2pu = NA,
      Krbc2pu = NA,
      kgutabs = NA,
      Kkidney2pu = NA,
      Klung2pu = NA,
      km = NA,
      Kmuc2air = NA,
      Kliver2pu = NA,
      Krest2pu = NA,
      Kblood2air = NA,
      kUrtc = NA,
      liver.density = NA,
      logHenry = NA,
      million.cells.per.gliver = NA,
      MW = NA,
      Pow = NA,
      pKa_Donor = NA,
      pKa_Accept = NA,
      MA = NA,
      Qcardiacc = NA,
      Qgfrc = NA,
      Qgutf = NA,
      Qliverf = NA,
      Qalvc = NA,
      Qkidneyf = NA,
      Qlungf = NA,
      Rblood2plasma = NA,
      Vgutc = NA,
      Vliverc = NA,
      Vartc = NA,
      Vkidneyc = NA,
      Vlungc = NA,
      vmax = NA,
      Vmucc = NA,
      Vvenc = NA,
      Vrestc = NA,
      KFsummary = NA,
      Fprotein.plasma = NA,
      fabs.oral = NA,
      Qgut_ = NA,
      Qintesttransport = NA
    ),
    param_units = alist(
      BW = "kg",
      Caco2.Pab = "1E-6 cm/s",
      Caco2.Pab.dist = "1E-6 cm/s",
      Clint = "uL/min/10^6 hepatocytes",
      Clint.dist = "uL/min/10^6 hepatocytes",
      Clmetabolismc = "L/h/kg BW",
      Funbound.plasma = "unitless fraction",
      Funbound.plasma.dist = "unitless fraction",
      Funbound.plasma.adjustment = "unitless coefficient",
      Fabsgut = "fraction",
      Fhep.assay.correction = "fraction",
      hematocrit = "percent volume RBCs in blood",
      Kgut2pu = "unitless ratio",
      Krbc2pu = "unitless ratio",
      kgutabs = "rate (1/hr)",
      Kkidney2pu = "unitless ratio",
      Klung2pu = "unitless ratio",
      km = "Michaelis-Menten concentration of half-maximal activity",
      Kmuc2air = "unitless ratio",
      Kliver2pu = "unitless ratio",
      Krest2pu = "unitless ratio",
      Kblood2air = "unitless ratio",
      kUrtc = "L/h/kg BW^(3/4)",
      liver.density = "g/cm^3",
      logHenry = "log10(atmosphers*m^3/mole)",
      million.cells.per.gliver = "cells/g liver",
      MW = "g/mol",
      Pow = "octanol:water partition coefficinet",
      pKa_Donor = "logarithmic",
      pKa_Accept = "logarithmic",
      MA = "phospholipid:water distribution coefficient",
      Qcardiacc = "L/h/kg BW^(3/4)",
      Qgfrc = "fraction",
      Qgutf = "fraction",
      Qliverf = "fraction",
      Qalvc = "L/h/kg BW^(3/4)",
      Qkidneyf = "fraction",
      Qlungf = "fraction",
      Rblood2plasma = "unitless ratio",
      Vgutc = "L/kg BW",
      Vliverc = "L/kg BW",
      Vartc = "L/kg BW",
      Vkidneyc = "L/kg BW",
      Vlungc = "L/kg BW",
      vmax = "max reaction velocity 1/min",
      Vmucc = "L/kg BW",
      Vvenc = "L/kg BW",
      Vrestc = "L/kg BW",
      KFsummary = "unitless",
      Fprotein.plasma = "fraction",
      fabs.oral = "fraction",
      Qgut_ = "fraction",
      Qintesttransport = "fraction"
    ),
    restrictive = TRUE) {

  stopifnot(c("Chemical", "Species") %in% names(data))

  this_species = data[["Species"]]
  this_chemical = data[["Chemical"]]

  # parameter names
  param_name <- c(
    "BW",
    "Caco2.Pab",
    "Caco2.Pab.dist",
    "Clint",
    "Clint.dist",
    "Clmetabolismc",
    "Funbound.plasma",
    "Funbound.plasma.dist",
    "Funbound.plasma.adjustment",
    "Fabsgut",
    "Fhep.assay.correction",
    "hematocrit",
    "Kgut2pu",
    "Krbc2pu",
    "kgutabs",
    "Kkidney2pu",
    "Klung2pu",
    "km",
    "Kmuc2air",
    "Kliver2pu",
    "Krest2pu",
    "Kblood2air",
    "kUrtc",
    "liver.density",
    "logHenry",
    "million.cells.per.gliver",
    "MW",
    "Pow",
    "pKa_Donor",
    "pKa_Accept",
    "MA",
    "Qcardiacc",
    "Qgfrc",
    "Qgutf",
    "Qliverf",
    "Qalvc",
    "Qkidneyf",
    "Qlungf",
    "Rblood2plasma",
    "Vgutc",
    "Vliverc",
    "Vartc",
    "Vkidneyc",
    "Vlungc",
    "vmax",
    "Vmucc",
    "Vvenc",
    "Vrestc",
    "KFsummary",
    "Fprotein.plasma",
    "fabs.oral",
    "Qgut_",
    "Qintesttransport"
  )

  # Default lower and upper bounds
  lower_bound_default <- alist(
    Clint = 0,
    Funbound.plasma = 0
  )
  upper_bound_default <- alist(
    Clint = 1000,
    Funbound.plasma = 1
  )

  # Fill missing lower and upper bounds
  lower_bound_missing <- base::setdiff(names(lower_bound_default), names(lower_bound))
  upper_bound_missing <- base::setdiff(names(upper_bound_default), names(upper_bound))

  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]

  # initialize optimization and use
  optimize_param <- rep(FALSE, length(param_name))
  optimize_param[param_name %in% c("Clint", "Funbound.plasma")] <- TRUE
  use_param <- rep(TRUE, length(param_name))

  param_units_vect <- sapply(param_units,
                             rlang::eval_tidy,
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  param_units_vect <- param_units_vect[param_name]

  lower_bound_vect <- sapply(lower_bound,
                             rlang::eval_tidy,
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  lower_bound_vect <- lower_bound_vect[param_name]

  upper_bound_vect <- sapply(upper_bound,
                             rlang::eval_tidy,
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  upper_bound_vect <- upper_bound_vect[param_name]

  par_DF <- data.frame("param_name" = param_name,
                       "param_units" = param_units_vect,
                       "optimize_param" = optimize_param,
                       "use_param" = use_param,
                       "lower_bound" = lower_bound_vect,
                       "upper_bound" = upper_bound_vect)

  # now get starting values
  par_DF <- get_starts_httk_gas_pbtk(data = data,
                                     par_DF = par_DF,
                                     this_species = this_species,
                                     this_chemical = this_chemical,
                                     restrictive = restrictive)

  return(par_DF)
}
