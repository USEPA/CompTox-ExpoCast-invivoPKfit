#' Calculates AUC for `httk`'s `gas_pbtk` PBPK model
#'
#' Calculated plasma concentration AUC vs time according to the `gas_pbtk`
#'
#' @inheritSection cp_httk_gas_pbtk Required parameters
#'
#' @inheritParams cp_httk_gas_pbtk
#'
#' @return A vector of blood or plasma AUC values  corresponding
#'  to `time`.
#'
#' @author Gilberto Padilla Mercado
#' @export auc_httk_gas_pbtk
#' @family built-in model functions
#' @family httk model functions
#' @family model auc functions
#'
auc_httk_gas_pbtk <- function(params, time, dose, route, medium = "plasma",
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
                    doses.per.day = NULL, # single dose
                    times = these_times,
                    iv.dose = this_route %in% "iv",
                    restrictive.clearance = restrictive,
                    input.units = "mg/kg",
                    output.units = "mg/L",
                    atol = 1E-9,
                    maxiter = 1E5,
                    small.time = 1E-6,
                    exp.conc = 0,
                    recalc.blood2plasma = FALSE,
                    recalc.clearance = FALSE,
                    suppress.messages = TRUE
                  )
          ) |>
            as.data.frame() |>
            subset(select = c("time", "AUC"))

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
        Time = x$Time,
        AUC = tmp[matched_times, "AUC"]
      )

      tmp$Route <- this_route
      tmp$Dose <- this_dose
      tmp$Medium <- this_medium

      tmp
    }
  )

  out_premerge <- do.call("rbind", res)
  rownames(out_premerge) <- NULL

  out <- merge(full_df, out_premerge)

  medium <- out$Medium
  out <- out[["AUC"]]

  # Transform data if Medium == "blood"
  out <- ifelse(medium %in% "blood", params$Rblood2plasma * out, out)

  return(out)
}
