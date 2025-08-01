#' Calculates plasma concentration for `httk`'s `gas_pbtk` model
#'
#' Calculated plasma concentrations vs time according to the `gas_pbtk` httk model
#'
#' @section Required parameters:
#' These are given by [httk::parameterize_gas_pbtk()].
#' Furthermore, they are transformed to a vector during the prefitting process.
#' The optimized parameters are `Clint` and `Funbound.plasma`. Because
#' these optimized parameters impact `Clmetabolismc`, `Krbc2pu`, `Rblood2plasma`
#' and `Fabsgut`, these are recalculated at the beginning of this function.
#'
#'
#' @param params A named numeric vector of model parameter values.
#' @param time A numeric vector of times, reflecting the time point when
#'  concentration is measured after the corresponding single bolus dose. Must be
#'  same length as `dose` and `iv.dose`, or length 1.
#' @param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `iv.dose`, or
#'  length 1. In this model, it is expected that this value represents a measurement
#'  of radioactive particles from a radiolabeling experiment.
#' @param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#' @param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as `time` and `dose`, or length 1.
#' @param this_chem A character vector naming the chemical for calculations in `httk`.
#' @param this_species A character vector naming the species for calculations in `httk`.
#' @param restrictive A logical value (TRUE or FALSE. Default: FALSE) that says whether the
#' assumption is that the clearance is restrictive or non-restrictive
#'
#' @return A vector of blood or plasma concentration values  corresponding
#'  to `time`.
#'
#' @author Gilberto Padilla Mercado
#' @export cp_httk_gas_pbtk
#' @family built-in model functions
#' @family httk model functions
#' @family model concentration functions
#'
cp_httk_gas_pbtk <- function(params, time, dose, route, medium = "plasma",
                             this_chem = NULL, this_species = NULL,
                             restrictive = TRUE) {
  # Make params into a list format
  params <- recalculate_httk_pbtk_params(params)
  # Create a data.frame to "track" the times
  full_df <- data.frame(
    Time = time/24, # convert to days
    Dose = dose,
    Route = route,
    Medium = medium
  )

  # Create another data.frame to get UNIQUE times
  uniq_df <- dplyr::distinct(full_df)
  rownames(uniq_df) <- NULL # Reset the rownames, with merge later it won't matter

  # Ensure this_chem and this_species is a single element vector
  if (length(this_chem) > 1) this_chem <- unique(this_chem)
  if (length(this_species) > 1) this_species <- unique(this_species)

  # Assert that it must be length 1 and not NULL
  stopifnot(!is.null(this_chem), !is.null(this_species),
            length(this_chem) == 1, length(this_chem) == 1)

  # Split-lapply pattern
  group_fct <- factor(with(uniq_df, paste(Dose, Route, Medium)))
  uniq_list <- split(uniq_df, group_fct) |> unname()

  # Solve httk function using these updated parameters
  res <- lapply(
    uniq_list,
    \(x) {
      this_dose <- unique(x$Dose)
      this_route <- unique(x$Route)
      this_medium <- unique(x$Medium)
      these_times <- unique(x$Time)

      tmp <- tryCatch(
        expr = {
          do.call(what = httk::solve_gas_pbtk,
                  args = list(
                    parameters = params,
                    species = this_species,
                    dose = this_dose,
                    times = these_times,
                    iv.dose = this_route %in% "iv",
                    restrictive.clearance = restrictive,
                    input.units = "mg/kg",
                    output.units = "mg/L",
                    default.to.human = FALSE,
                    exp.conc = 0,
                    doses.per.day = NULL, # single dose
                    recalc.blood2plasma = FALSE,
                    recalc.clearance = FALSE,
                    suppress.messages = TRUE
                  )
          ) |>
            as.data.frame() |>
            subset(select = c("time", "Cplasma"))

        }, warning = function(w) {
          message("Here is a warning")
          len_time <- length(unique(c(0, x$Time)))
          data.frame(time = unique(c(0, x$Time)),
                     Cplasma = rep(-1, len_time))

        }, error = function(e) {

          message("Here is an error")
          # if there is an error, return negative values
          # this ensures this parameter set is discarded during optimization
          len_time <- length(unique(c(0, x$Time)))
          data.frame(time = unique(c(0, x$Time)),
                     Cplasma = rep(-1, len_time))
        }
      )

      matched_times <- sapply(
        x$Time,
        \(z) {
          which.min(abs(z - tmp$time))
        }
      )

      tmp <- data.frame(
        Time = tmp[matched_times, "time"],
        Cplasma = tmp[matched_times, "Cplasma"]
      )

      tmp$Route <- this_route
      tmp$Dose <- this_dose
      tmp$Medium <- this_medium

      tmp
    }
  )

  out_premerge <- do.call("rbind", res)
  rownames(out_premerge) <- NULL

  out <- dplyr::left_join(full_df, out_premerge, by = c("Time", "Dose", "Route", "Medium"))

  medium <- out$Medium
  out <- out[["Cplasma"]]

  # Transform data if Medium == "blood"
  out <- ifelse(medium %in% "blood", params$Rblood2plasma * out, out)

  return(out)
}
