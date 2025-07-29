#' Recalculating parameters for `httk`'s pbtk or gas_pbtk model
#'
#' This set of functions recalculates the parameters that change with each new
#' set of Funbound.plasma and Clint values. In summary these are:
#' metabolic clearance (Clmetabolismc), red blood cell partitioning coefficient
#' (krbc2pu), blood-to-plasma ratio (Rblood2plasma), and fraction absorbed by
#' the gut (Fabsgut).
#'
#' @section Clmetabolismc:
#' The formula to recalculate this parameters is as follows:
#' \deqn{Clmetabolismc = Clint\times million.cells.per.gliver\times Vliverc\times
#' \frac{10^3g}{kg}\times liver.density\times Funbound.plasma\times
#' \frac{60min}{hr}\times \frac{L}{10^6\mu L}}
#'
#' Note that because Clmetabolismc is in \eqn{\frac{\mu L}{min\cdot 10^6 cells}},
#' Vliverc is in \eqn{\frac{L}{kg\ BW}} and liver.density is in \eqn{\frac{kg}{L}},
#' there are some unit conversions included to give final units in
#' \eqn{\frac{L}{hr\cdot kg\ BW}}
#'
#'
#' @section krbc2pu:
#' The formula to recalculate this parameter is as follows:
#' \deqn{krbc2pu = Fint\times Kint\times KAPPAcell2pu\times Fcell\times Kcell}
#' where Fint, Fcell, and Kcell are taken or calculated from values in [httk::tissue.data],
#' and KAPPAcell2pu is estimated once during [httk::calc_ionization()] and
#' [httk::parameterize_schmitt()].
#' \deqn{Kint = (1 - Fprotein.plasma + \left(\frac{0.37}{Funbound.plasma - (1-Fprotein.plasma)}\right)}
#' Because only Kint needs to be recalculated, the rest is saved as a summary
#' constant called KFsummary and the other additional parameter included is Fprotein.plasma.
#'
#'
#' @section Rblood2plasma:
#' The formula to calculate this parameter is:
#' \deqn{Rblood2plasma = 1 - hematocrit + (hematocrit\times krbc2pu\times Funbound.plasma)}
#'
#'
#' @section Fabsgut:
#' The formula to calculate this parameter is:
#' \deqn{Fabsgut = fabs.oral\times fgut.oral}
#' where
#' \deqn{fabs.oral = \mathrm{min}(1, 1 - \left(1 + 2\times peff \times \frac{MRT\times Rsi}{7}
#' \times \frac{60}{10^4}\right)^{-7})}
#'
#' \deqn{fgut.oral = \mathrm{min}(1,
#' \frac{Qgut}{Qgut + \left(Funbound.plasma\times \frac{Clmetabolismc\times BW}{100}\right)}\times
#' \frac{Qgut}{Qintesttransport + Qgut})}
#' and
#' \deqn{Qgut = \frac{Qvilli\times CLperm}{Qvilli + CLperm}}
#' \deqn{Qvilli = Qvillif\times Qsif\times Qgutf\times Qcardiacc\times BW^{3/4}}
#' \deqn{CLperm = p_{eff}\times Asi\times \left(\frac{1000}{10^4}\times 100\times 3600\right)}
#' \deqn{p_{eff} = 10^{0.4926\times \mathrm{log10}(Caco2.Pab) - 0.1454} \text{ and }
#' \frac{p_{eff} + 1.815}{1.039}\text{ when Species == "rat"}}
#' \deqn{Asi = 0.66\times BW/70 \text{ or } 71/100^2 \text{ when Species == "rat"}}
#' \deqn{Qintesttransport = 0.1\times\left(\frac{BW}{70}\right)^{3/4}}
#'
#'
#' @param params A list of parameter = value pairs that will be used in calculations.
#' @param dtxsid The DTXSID of a chemical, by default taken from the 'Chemical' column in the data.
#' @param species The species of a subject, by default taken from the 'Species' column in the data.
#'
#' @returns An updated list of parameters
#' @export
#'
recalculate_httk_pbtk_params <- function(params, dtxsid, species) {

  KFsummary <- Fprotein.plasma <- fabs.oral <- hematocrit <-  Qgut_ <- Qintesttransport <- NULL
  Clint <- million.cells.per.gliver <- Vliverc <- liver.density <- Funbound.plasma <- BW <- NULL

  list2env(as.list(params), envir = as.environment(-1))

  Clmetabolismc <- Clint * million.cells.per.gliver *
    Vliverc * 1E3 * liver.density * Funbound.plasma * 60 * 1E-6

  Kint <- 1 - Fprotein.plasma + (0.37 / Funbound.plasma - (1 - Fprotein.plasma))
  Krbc2pu <- KFsummary * Kint

  Rblood2plasma <- 1 - hematocrit + (hematocrit * Krbc2pu * Funbound.plasma)

  fgut.oral <- min(
    1,
    (Qgut_ / (Qgut_ + (Funbound.plasma * Clmetabolismc * BW / 100))) *
      (Qgut_ / (Qintesttransport + Qgut_))
  )

  Fabsgut <- fabs.oral * fgut.oral

  params[["Clmetabolismc"]] <- Clmetabolismc
  params[["Krbc2pu"]] <- Krbc2pu
  params[["Rblood2plasma"]] <- Rblood2plasma
  params[["Fabsgut"]] <- Fabsgut
  if (!nzchar(params[["Caco2.Pab.dist"]])) params[["Caco2.Pab.dist"]] <- NA_character_
  if (!nzchar(params[["Clint.dist"]])) params[["Clint.dist"]] <- NA_character_
  if (!nzchar(params[["Funbound.plasma.dist"]])) params[["Funbound.plasma.dist"]] <- NA_character_

  # For whatever reason, additional parameters trip up the solver
  params[c("KFsummary", "Fprotein.plasma",
           "fabs.oral", "Qgut_", "Qintesttransport")] <- NULL

  # Similar to httk:::set_httk_precision
  params <- lapply(params, \(x) ifelse(is.numeric(x), round(signif(x, 4), 9), x))

  return(params)

}
