#' Fill parameters for 1-compartment model
#'
#' This is used to "fill" parameters prior to running `httk::solve_3compartment2()`.
#' During pre-fitting, the pKa_Donor & pKa_Accept parameters may hold multiple
#' values and to keep data tidy, each value was split into it's own row. This
#' function reconvenest these values into the appropriate format and outputs a
#' list which `httk` can use.
#'
#' @param params Named numeric vector of parameters for the 1-compartment model
#' @return A named list of parameters,to be used with `httk`'s
#'  `3compartment2` model.
#' @author Gilberto Padilla Mercado

fill_params_httk_3comp2 <- function(params) {

  # Make params that are not pKa a list
  not_pKa <- params[-grep("pKa_(Accept|Donor)_\\d", names(params))] |> as.list()

  # Reconvene pKa params (grouped by Accept/Donor)
  par_pKa_Accept <- params[grep("pKa_Accept_\\d", names(params))]
  par_pKa_Donor <- params[grep("pKa_Donor_\\d", names(params))]

  par_pKa_Accept <- list(
    pKa_Accept = c(pKa_Accept = num2string(par_pKa_Accept))
  )

  par_pKa_Donor <- list(
    pKa_Donor = c(pKa_Donor = num2string(par_pKa_Donor))
  )

  # Concatenate all parameters into one list
  pars_out <- c(not_pKa, par_pKa_Accept, par_pKa_Donor)

  return(pars_out)
}
