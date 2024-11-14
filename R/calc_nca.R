#' Non-compartmental analysis
#'
#' Do non-compartmental analysis on a single-dose set of concentration vs. time
#' data
#'
#' This function is a wrapper around [PK::nca()] to do non-compartmental
#' analysis, after automatically detecting the study design. It additionally
#' calls [get_peak()] to calculate the peak concentration and time of peak
#' concentration.
#'
#' # Automatic detection of study design
#'
#' [PK::nca()] understands three different study designs, and requires the user
#' to specify which one is being used.
#'
#' - `ssd`: Serial sampling design. Each observation is from a different subject.
#' - `complete`: Every subject was observed at every time point.
#' - `batch`: Each subject was observed at multiple time points, but not at every time point.
#'
#' To automatically detect which study design is applicable, this function first
#' sorts the data by increasing time. Then, a table of time vs. series ID is
#' created, with 1 indicating that a measurement exists for the corresponding
#' time point/series ID combination, and 0 indicating that a measurement does
#' not exist. If the column sums of this table are all 1, then it is a serial
#' sampling design, except if there is only one observation per time point, it
#' is a complete sampling design, and if there are multiple observations for
#' some time points and only one observation for other time points, it is a
#' batch design. If the column sums are all equal to the number of rows of the
#' table, then it is a complete sampling design. Otherwise, it is a batch
#' sampling design.
#'
#' # Parameters estimated by NCA
#'
#' - `AUC_infinity`: The area under the concentration-time curve, extrapolated out to infinite time. Estimated using the trapezoidal rule, with a tail area correction calculated using the slope of the last 3 data points (by default).
#' - `AUC_tlast`: The area under the concentration-time curve, calculated at the last observed time point. Estimated using the trapezoidal rule.
#' - `AUMC_infinity`: The area under the concentration-time first moment curve (the area under the AUC vs. time), extrapolated out to infinite time. Estimated using the trapezoidal rule, with a tail area correction calculated using the slope of the last 3 data points (by default).
#' - `CLtot`: The total clearance rate.  Only calculated for `route == 'iv'`. If `route == 'oral'`, this is `NA_real_`, and only `CLtot/Fgutabs` is calculated.
#' - `CLtot/Fgutabs`: The total clearance rate, normalized by the oral bioavailability. Only calculated for `route == 'oral'`. If `route == 'iv'`, this is `NA_real_`, and only `CLtot/` is calculated.
#' - `Cmax`: The peak concentration. For `route == 'iv'`, this is expected to be the concentration at the earliest time; for `route == 'oral'`, it is not. This and `tmax` are calculated using [get_peak()], not by [PK::nca()].
#' - `halflife`: The half-life of elimination.  Only calculated for `route == 'iv'`. If `route == 'oral'`, this is `NA_real_`, because half-life estimates are not valid for oral data.
#' - `MRT`: The mean residence time. Only calculated for `route == 'iv'`. If `route == 'oral'`, this is `NA_real_`, and only `MTT` is calculated.
#' - `MTT`: The mean transit time (the sum of MRT and mean absorption time). Only calculated for `route == 'oral'`. If `route == 'iv'`, this is `NA_real_`, and only `MRT` is calculated.
#'  -`tmax`: The time of peak concentration. For `route == 'iv'`, this is expected to be the earliest time; for `route == 'oral'`, it is not.  This and `Cmax` are calculated using [get_peak()], not by [PK::nca()].
#' - `Vss`: The volume of distribution at steady state (`AUMC_infinity/AUC_infinity^2`). If `route == 'oral'`, this is `NA_real_`, because `Vss` estimates are not valid for oral data.
#'
#' # Output
#'
#' The output is a data.frame with 9 rows (one for each NCA parameter) and a
#' number of variables equal to `length(method) + 3`.
#'
#' The variables are
#'
#'  - `design`: The automatically-detected design. One of `ssd`, `complete`, or `batch` (or `NA_character_` if no analysis could be done).
#'  - `param_name`: The name of each NCA parameter.
#'  - `param_value`: The value of each NCA parameter.
#'  - `param_sd_[method]`: The parameter standard error estimated by the corresponding method.
#'
#'
#'
#' @param time A numeric vector of time points.
#' @param conc A numeric vector of concentrations. If detected (above limit of
#'   detection/quantification), contains the measured value; if not detected
#'   (below LOD/LOQ), contains the LOD/LOQ.
#' @param detect A logical vector: Whether each concentration was detected
#'   (above LOD/LOQ) or not.
#' @param series_id Optional: A variable that can be coerced to a factor,
#'   identifying individual time series (e.g., individual replicates --
#'   individual subjects, or replicate dose groups). Default NULL, in which case
#'   each observation will be assumed to have a different series ID. In other
#'   words, a serial sampling design will be assumed, in which each observation
#'   is from a different subject.
#' @param dose A numeric scalar: The dose for this data set.
#' @param route A character scalar: The route of administration for this data
#'   set. Currently, only "oral" and "iv" are supported.
#' @param method As for [PK::nca()]: the method to use for calculation of
#'   confidence intervals. Default `'z'` (this differs from the [PK::nca()]
#'   default).
#' @param ... Other arguments that will be passed to [PK::nca()] (other than
#'   `data`, `design`, and `method`: *i.e.*, `n.tail`, `nsample`)
#' @return A `data.frame` with 9 rows and `length(method) + 3` variables. See
#'   Details.
#' @export
#' @import PK
#' @author Caroline Ring
calc_nca <- function(time,
                    conc,
                    detect,
                   series_id = NULL,
                   dose,
                   route,
                   method = "z",
                   ...) {

  if (length(time) > 0 && length(conc) > 0 && length(dose) > 0 &&
      anyNA(c(time, conc, dose)) && anyNA(detect)) {

    dose <- unique(dose)

    if (is.null(series_id)) {
      series_id <- rep(NA_integer_, length(conc))
    }

    # order everything by increasing time
    ord <- order(time)
    time <- time[ord]
    conc <- conc[ord]
    detect <- detect[ord]
    series_id <- series_id[ord]

    # Calculate area under the concentration-time curve for NCA.
    # At present we do not know individual animal IDs or animal-group IDs for each point,
    # so we will have to just assume that each time point is a different animal.
    # substitute any nondetects with LOQ/2. This is not sophisticated but it is fast.

    conc <- ifelse(detect %in% FALSE, conc * 0.5, conc)

    if (all(is.na(series_id))) {
      # assume every obs is a different animal
      series_id <- seq_along(conc)
    }


    ntab <- table(time, series_id)
    m <- matrix(ntab, ncol = ncol(ntab))

    obs_per_time <- table(time)
    obs_per_seriesID <- table(series_id)
    times_per_seriesID <- colSums(table(time, series_id))

    # if there is only one observation per time
    # then pretend all observations are from the same id
    if (all(obs_per_time %in% 1)) {
      series_id <- rep(1L, length(time))
    }

    obs_per_time <- table(time)
    obs_per_seriesID <- table(series_id)
    times_per_seriesID <- colSums(table(time, series_id))

    ntab <- table(time, series_id)
    m <- matrix(ntab, ncol = ncol(ntab))

    if (all(m == 1)) {
      # every subject measured at every time point
      design <- "complete"
    } else if (any(obs_per_time == 1) && any(obs_per_time > 1)) {
      # if 1 observation for some time point and multiple observations for others,
      # it should be "batch" but it will fail due to a bug in PK::nca.batch()
      # munge time very slightly to make 1 observation for each time point,
      # and use an ssd design
      time_split <- split(time, time)
      min_time_nz <- min(time[time > 0])
      time_split <- sapply(time_split,
                           function(this_time) {
                             if (length(this_time) == 1) {
                               this_time
                             } else {
                               # add a small fuzz factor: 1% of smallest time
                               this_time + runif(length(this_time),
                                                 min = 0,
                                                 max = 0.01 * min_time_nz)
                             }
                           })
      time <- unsplit(time_split, time)
      design <- "complete"
    } else if (all(obs_per_time > 1) && length(unique(obs_per_seriesID)) > 1) {
      # multiple observations for every time point, but not all subjects at every time point
      # note this will break if there is 1 time point for some subjects and multiple for others
      if (any(times_per_seriesID < 2)) {
        # pretend every obs is a different series id
        series_id <- seq_along(conc)
        design <- "ssd"
      } else {
        # if there is more than one time point measured for each subject,
        # and more than one obs per time point,
        # but not all time points measured for all subjects
        design <- "batch"
      }
    } else if (length(unique(obs_per_seriesID)) == 1) {
      # one obs per subject per time point
      design <- "ssd"
    } else {
      design <- "batch"
    }

    data <- data.frame(id = series_id, conc = conc, time = time)

    pk_out <- tryCatch(
      {
        tmp <- suppressMessages(
          suppressWarnings(
            do.call(PK::nca,
                    args = c(list(data = data,
                                  dose = dose,
                                  design = design,
                                  method = method),
                             list(...)))
          )
        )
        tmp_est <- tmp$est[, 1]
        tmp_se <- tmp$CIs[, c("stderr", "method")]
        # if more than one method, reshape the output to have one stderr column per method
        tmp_se_list <- sapply(
          method,
          function(this_method) {
            this_tmp <- tmp_se[tmp_se[, "method"] %in% this_method, 1]
            this_tmp
          })
        tmp_out <- cbind(tmp_est, tmp_se_list)
        tmp_out
      },
      error = function(err) {
        tmp <- matrix(nrow = 7,
                      ncol = length(method) + 1)
        rownames(tmp) <- c("AUC_tlast",
                           "AUC_infinity",
                           "AUMC_infinity",
                           "MRT",
                           "halflife",
                           "CLtot",
                           "Vss")
        colnames(tmp) <- c("est",
                           paste("se",
                                 method,
                                 sep = "."))
        tmp
      })
  } else { # If data are zero
    design <- NA_character_
    pk_out <- matrix(nrow = 7,
                     ncol = length(method) + 1)
    rownames(pk_out) <- c("AUC_tlast",
                          "AUC_infinity",
                          "AUMC_infinity",
                          "MRT",
                          "halflife",
                          "CLtot",
                          "Vss")
    colnames(pk_out) <- c("est",
                          paste("se",
                                method,
                                sep = "."))
  }

  if (all(route %in% "oral")) {
    rownames(pk_out) <- c("AUC_tlast",
                          "AUC_infinity",
                          "AUMC_infinity",
                          "MTT",
                          "halflife",
                          "CLtot/Fgutabs",
                          "Vss")
    # halflife and Vss are not valid under oral administration, per ?PK::nca
    pk_out["halflife", ] <- NA_real_
    pk_out["Vss", ] <- NA_real_
    # fill in CLtot and MRT as NA
    pk_out <- rbind(pk_out,
                    "CLtot" = rep(NA_real_, ncol(pk_out)),
                    "MRT" = rep(NA_real_, ncol(pk_out)))
  } else {
    rownames(pk_out) <- c("AUC_tlast",
                          "AUC_infinity",
                          "AUMC_infinity",
                          "MRT",
                          "halflife",
                          "CLtot",
                          "Vss")
    # fill in oral-only params as NA
    pk_out <- rbind(pk_out,
                    "CLtot/Fgutabs" = rep(NA_real_, ncol(pk_out)),
                    "MTT" = rep(NA_real_, ncol(pk_out)))
  }


  # also compute tmax, Cmax
  peak <- unlist(get_peak(x = time, y = conc))

  pk_out <- rbind(pk_out,
                  "tmax" = c(peak[1], rep(NA_real_, ncol(pk_out) - 1)),
                  "Cmax" = c(peak[2], rep(NA_real_, ncol(pk_out) - 1))
  )

  # convert to data.frame
  outval <- as.data.frame(pk_out)
  names(outval) <- c("param_value", paste0("param_sd_", method))
  outval$param_name = rownames(pk_out)
  outval$design <- design

  # put columnsin right order
  outval <- outval[, c("design", "param_name", "param_value", paste0("param_sd_", method))]

  # ensure parameters are sorted alphabetically
  outval <- outval[order(outval$param_name), ]

  # remove rownames
  rownames(outval) <- NULL

  outval
}
