#' Get parameters to fit `httk`'s `3compartment2` PBPK model
#'
#' Get parameters to fit the `3compartment2`
#' model from the `httk` package (Wambaugh, Schacht, and Ring. 2025).
#'
#' @section Required parameters:
#' These are given by the \link[httk]{parameterize_3comp2} function in `httk`.
#' Furthermore, they are transformed to a vector during hte prefitting process.
#'
#' @param data The data set to be fitted (e.g. the result of [preprocess_data()])
#' @param lower_bound A mapping specified using a call to [ggplot2::aes()],
#'  giving the lower bounds for each variable, as expressions which may include
#'  variables in `data`.
#' @param upper_bound A mapping specified using a call to [ggplot2::aes()],
#'  giving the upper bounds for each variable, as expressions which may include
#'  variables in `data`.
#' @param param_units A mapping specified using a call to [ggplot2::aes()],
#'  giving the units for each variable, as expressions which may include
#'  variables in `data`.
#' @param this_chemical A character vector naming the chemical for calculations in `httk`.
#' @param this_species A character vector naming the species for calculations in `httk`.
#' @param restrictive A logical value (TRUE or FALSE. Default: FALSE) that says whether the
#' assumption is that the clearance is restrictive or non-restrictive
#'
#' @return A vector of blood or plasma concentration values  corresponding
#'  to `time`.
#' @author Gilberto Padilla Mercado
#' @export get_params_httk_3comp2
#' @family built-in model functions
#' @family httk model functions
#' @family model concentration functions
#'
get_params_httk_3comp2 <- function(
    data,
    lower_bound = ggplot2::aes(
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
      Kliver2pu = NA,
      Krest2pu = NA,
      Kblood2air = NA,
      liver.density = NA,
      logHenry = NA,
      million.cells.per.gliver = NA,
      MW = NA,
      Pow = NA,
      # pKa_Donor and pKa_Accept are variable length but will be included in final
      MA = NA,
      Qcardiacc = NA,
      Qgfrc = NA,
      Qgutf = NA,
      Qliverf = NA,
      Qalvc = NA,
      Rblood2plasma = NA,
      Vgutc = NA,
      Vliverc = NA,
      Vrestc = NA
    ),
    upper_bound = ggplot2::aes(
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
      Kliver2pu = NA,
      Krest2pu = NA,
      Kblood2air = NA,
      liver.density = NA,
      logHenry = NA,
      million.cells.per.gliver = NA,
      MW = NA,
      Pow = NA,
      # pKa_Donor and pKa_Accept are variable length but will be included in final
      MA = NA,
      Qcardiacc = NA,
      Qgfrc = NA,
      Qgutf = NA,
      Qliverf = NA,
      Qalvc = NA,
      Rblood2plasma = NA,
      Vgutc = NA,
      Vliverc = NA,
      Vrestc = NA
    ),
    param_units = ggplot2::aes(
      BW = "kg",
      Caco2.Pab = "1E-6 cm/s",
      Caco2.Pab.dist_m = "1E-6 cm/s",
      Caco2.Pab.dist_u = "1E-6 cm/s",
      Caco2.Pab.dist_l = "1E-6 cm/s",
      Clint = "uL/min/10^6 hepatocytes",
      Clint.dist_m = "uL/min/10^6 hepatocytes",
      Clint.dist_l = "uL/min/10^6 hepatocytes",
      Clint.dist_u = "uL/min/10^6 hepatocytes",
      Clmetabolismc = "L/h/kg BW",
      Funbound.plasma = "unitless fraction",
      Funbound.plasma.dist_m = "unitless fraction",
      Funbound.plasma.dist_l = "unitless fraction",
      Funbound.plasma.dist_u = "unitless fraction",
      Funbound.plasma.adjustment = "unitless coefficient",
      Fabsgut = "fraction",
      Fhep.assay.correction = "fraction",
      hematocrit = "percent volume RBCs in blood",
      Kgut2pu = "unitless ratio",
      Krbc2pu = "unitless ratio",
      kgutabs = "unitless ratio",
      Kliver2pu = "unitless ratio",
      Krest2pu = "unitless ratio",
      Kblood2air = "unitless ratio",
      liver.density = "g/cm^3",
      logHenry = "log10(atmosphers*m^3/mole)",
      million.cells.per.gliver = "cells/g liver",
      MW = "g/mol",
      Pow = "octanol:water partition coefficinet",
      # pKa_Donor, pKa_Accept,
      MA = "phospholipid:water distribution coefficient",
      Qcardiacc = "L/h/kg BW^(3/4)",
      Qgfrc = "fraction",
      Qgutf = "fraction",
      Qliverf = "fraction",
      Qalvc = "fraction",
      Rblood2plasma = "unitless ratio",
      Vgutc = "L/kg BW",
      Vliverc = "L/kg BW",
      Vrestc = "L/kg BW"
    ),
    this_species,
    this_chemical,
    restrictive = TRUE) {

  # parameter names
  param_name <- c("BW",
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
                  "Kliver2pu",
                  "Krest2pu",
                  "Kblood2air",
                  "liver.density",
                  "logHenry",
                  "million.cells.per.gliver",
                  "MW",
                  "Pow",
                  # pKa values are added in get_starts
                  "MA",
                  "Qcardiacc",
                  "Qgfrc",
                  "Qgutf",
                  "Qliverf",
                  "Qalvc",
                  "Rblood2plasma",
                  "Vgutc",
                  "Vliverc",
                  "Vrestc"
  )

  # Default lower and upper bounds
  lower_bound_default <- ggplot2::aes(
    Clint = 0,
    Funbound.plasma = 0
  )
  upper_bound_default <- ggplot2::aes(
    Clint = 500,
    Funbound.plasma = 1
  )

  # Fill missing lower and upper bounds
  lower_bound_missing <- setdiff(names(lower_bound_default), names(lower_bound))
  upper_bound_missing <- setdiff(names(upper_bound_default), names(upper_bound))

  lower_bound[lower_bound_missing] <- lower_bound_default[lower_bound_missing]
  upper_bound[upper_bound_missing] <- upper_bound_default[upper_bound_missing]

  # initialize optimization and use
  optimize_param <- rep(FALSE, length(param_name))
  optimize_param[param_name %in% c("Clint", "Funbound.plasma")] <- TRUE
  use_param <- rep(TRUE, length(param_name))

  param_units_vect <- sapply(param_units,
                             \(x) { rlang::eval_tidy(x, data = data) },
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  param_units_vect <- param_units_vect[param_name]

  lower_bound_vect <- sapply(lower_bound,
                             \(x) { rlang::eval_tidy(x, data = data) },
                             simplify = TRUE,
                             USE.NAMES = TRUE)
  lower_bound_vect <- lower_bound_vect[param_name]

  upper_bound_vect <- sapply(upper_bound,
                             \(x) { rlang::eval_tidy(x, data = data) },
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
  par_DF <- get_starts_httk_3comp2(data = data,
                                   par_DF = par_DF,
                                   this_species = this_species,
                                   this_chemical = this_chemical,
                                   restrictive = restrictive) |>
    suppressWarnings()

  # check to ensure starting values are within bounds
  # if not, then replace them by a value halfway between bounds
  start_low <- (par_DF$start < par_DF$lower_bound) %in% TRUE
  start_high <- (par_DF$start > par_DF$upper_bound) %in% TRUE
  start_nonfin <- !is.finite(par_DF$start)

  par_DF[start_low | start_high | start_nonfin,
         "start"] <- rowMeans(
           cbind(par_DF[start_low | start_high | start_nonfin,
                        c("lower_bound",
                          "upper_bound")]
           )
         )

  return(par_DF)
}
