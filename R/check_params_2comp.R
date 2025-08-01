#' Check 2-compartment model parameters
#'
#' Check to make sure required parameters are present to evaluate 2-compartment
#' model for a given route and medium
#'
#' @param params A named numeric vector of parameters for the 2-compartment model.
#' @param route A character vector of routes: "iv" and/or "oral".
#' @param medium A character vector of tissue media: "plasma" and/or "blood".
#' @param ... Additional arguments (not currently used)
#' @return Character: A message. If all required parameters are present for the
#'   given media & routes, the message is "Parameters OK". If required
#'   parameters for the oral route are missing, the message is "Error: For
#'   2-compartment oral model, missing parameters (comma-separated list of
#'   parameter names)".  If required
#'   parameters for the IV route are missing, the message is "Error: For
#'   2-compartment oral model, missing parameters (comma-separated list of
#'   parameter names)".
#' @author Caroline Ring

check_params_2comp <- function(params,
                               route,
                               medium,
                               ...) {

  msg <- "Parameters OK"
  params <- unlist(params)

  # check for any missing parameters
  # required params for oral dose
  if (any(route %in% "oral")) {
    missing_params <- base::setdiff(c("kelim", "k21", "k12", "Fgutabs_V1", "kgutabs"),
                              names(params)[is.finite(params)])
    if (length(missing_params) > 0) {
      msg <- cli::cli_fmt(
        cli::cli_text("For 2-compartment oral model, missing parameters: {missing_params}")
      )
    }
  }

  # required params for IV dose
  if (any(route %in% "iv")) {
    missing_params <- base::setdiff(c("kelim", "V1", "k12", "k21"),
                              names(params)[is.finite(params)])
    if (length(missing_params) > 0) {
      msg <- cli::cli_fmt(
        cli::cli_text("For 2-compartment IV model, missing parameters: {missing_params}")
      )
    }
  }

  if (any(medium %in% "blood") && !"Rblood2plasma" %in% names(params)[is.finite(params)]) {
    msg <- cli::cli_fmt(
      cli::cli_text("For 2-compartment model with blood data, missing parameter: Rblood2plasma")
    )
  }

  return(msg)
}
