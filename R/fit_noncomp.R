fit_noncomp <- function(data) {

  ### PK package doesn't seem to like replicate timepoints
  # here is a function to to take the average concentrations of replicate timepoints
  calc_avg_conc <- function(data) {

    ### convert to data.frame, but fix this later because converting back and forth is silly
    data <- as.data.frame(data)

    ### if there are replicate timepoints, take the average Value for each timepoint
    if(length(unique(data$time)) < length(data$time)) {
      data <- data %>%
        dplyr::group_by(time) %>%
        dplyr::do(dplyr::mutate(., mean = mean(.$conc, na.rm = TRUE)))

      ### if all values are NA, output will return 'NaN's
      ### coerces 'NaN's back to NA
      data$mean[is.nan(data$mean)] <- NA

      ### remove replicate timepoints
      data <- data %>%
        dplyr::distinct(time, .keep_all = TRUE)

      ### delete Value column and reassign new mean column as Value column
      data$conc <- NULL
      names(data)[names(data) == "mean"] <- "conc"
    }

    data <- as.data.table(data)

    data
  }

  ### apply that function to unique combos
  data_list <- list()
  for(this_compound in unique(data[, compound])) {
    this_compound_data <- data[compound == this_compound]

    ### list of columns by which to split data.set
    col_list <- list(this_compound_data$reference,
                     this_compound_data$dose,
                     this_compound_data$route,
                     this_compound_data$media)

    ### split data.set into list of data.tables
    split_list <- split(this_compound_data, col_list)

    ### take average of replicate timepoints for each data.table
    split_list_avg <- lapply(split_list, calc_avg_conc)

    sub_data <- do.call(rbind, split_list_avg)

    data_list[[this_compound]] <- sub_data
  }

  ### combine data.tables back into one large data.set
  data <- do.call(rbind, data_list)

  ### normalize all concentrations by dose to fit jointly
  data[, conc_norm := conc / dose]

  pk_fit_table <- NULL

  for (this_cas in sort(unique(data$cas))) {
    this_subset <- subset(data, cas == this_cas & !is.na(conc_norm))

    for (this_species in sort(unique(this_subset$species))){

      # Need oral and iv to estimate bioavailability:
      if ("po" %in% unique(this_subset$route) & "iv" %in% unique(this_subset$route)) {

        this_iv_time_list <- list()
        this_iv_conc_list <- list()
        this_iv_id_list <- list()
        this_oral_time_list <- list()
        this_oral_conc_list <- list()
        this_oral_id_list <- list()
        this_row <- data.frame(compound = this_subset$compound[1],
                               cas = this_cas,
                               reference = NA,
                               # Source = NA,
                               species = this_species,
                               mean_iv_dose = mean(subset(this_subset, route == "iv")$dose), ### averaged doses for iv and po
                               mean_po_dose = mean(subset(this_subset, route == "po")$dose),
                               model_type = "noncompartmental",
                               param_value_type = "estimated",
                               auc_po = NA,
                               vd_po = NA,
                               cl_po = NA,
                               auc_iv = NA,
                               vd_iv = NA,
                               cl_iv = NA,
                               fbio = NA,
                               stringsAsFactors = F)

        this_row <- rbind(this_row, this_row)

        rownames(this_row) <- NULL

        this_row[2,"param_value_type"] <- "estimated std dev"

        for (this_route in unique(this_subset$route)) {

          this_route_subset <- subset(this_subset, route == this_route)

          # new.route.subset <- NULL
          # print(this.route.subset$Compound)
          this_route_subset$subject <- as.character(this_route_subset$subject) ### added because Subject was logical and characters were not coerced correctly

          if (any(is.na(this_route_subset$subject))) {

            for (this_reference in unique(this_route_subset$reference)) {

              for (this_dose in unique(this_route_subset$dose)) {
                this_route_subset[is.na(this_route_subset$subject) & this_route_subset$dose == this_dose, "subject"] <- paste(this_reference, this_route, this_dose, sep = "-")
              }
            }
          }

          this_route_subset$id <- paste(this_row[1, "reference"], this_route_subset$subject, sep = "-")

          ### conc and time are required fields for the PK package
          this_route_subset$conc <- this_route_subset$conc_norm
          this_route_subset$time <- this_route_subset$time

          for (this_id in unique(this_route_subset$id)) {

            ### comment here
            if (dim(subset(this_route_subset, id == this_id))[1] < 3) {

              this_route_subset <- subset(this_route_subset, id != this_id)
            }
          }

          ### order this.route.subset by time
          this_route_subset <- this_route_subset[order(this_route_subset$time), ]

          if (nrow(this_route_subset) > 0) {

            this_id_list <- list()

            for (this_id in unique(this_route_subset[, id])) {

              this_id_data <- this_route_subset[id == this_id]

              this_id_data <- as.data.frame(this_id_data)

              if (length(unique(this_id_data$time)) > 1) {

                this_id_data_group <- this_id_data %>%
                  dplyr::group_by(time) %>%
                  dplyr::tally()
                this_id_data <- dplyr::full_join(this_id_data, this_id_data_group)
                this_id_data <- as.data.table(this_id_data)
                this_id_data <- this_id_data[n == max(n)]
                this_id_data <- this_id_data[, n := NULL]
              }

              this_id_list[[this_id]] <- this_id_data
            }

            this_route_subset <- do.call(rbind, this_id_list)
          }

          if (dim(this_route_subset)[1] > 2) {

            if (this_route == "po") {

              dose_arg <- 0
            } else if (this_route == "iv") {

              dose_arg <- 1
            }

            ### depending on length time, try .complete or .batch
            if (length(this_route_subset$time) == length(unique(this_route_subset$time))) {

              out <- try(PK::nca.complete(data = this_route_subset, dose = dose_arg, method = "z"))
            } else {
              out <- try(PK::nca.batch(data = this_route_subset, dose = dose_arg, method = "z"))
            }

browser()
            if (class(out) != "try-error") {

              out <- out$cis

              if (!is.na(out[1, "est"] > 0)) if (out[1, "est"] > 0) {

                this_row[1, paste("auc", this_route, sep = ".")] <- out[1, "est"]
                this_row[2, paste("auc", this_route, sep = ".")] <- out[1, "stderr"]
              }

              if (!is.na(out[6, "est"] > 0)) if (out[6, "est"] > 0) {

                this_row[1, paste("cl", this_route, sep = ".")] <- out[6, "est"]
                this_row[2, paste("cl", this_route, sep = ".")] <- out[6, "stderr"]
              }

              if(!is.na(out[7, "est"] > 0)) if (out[7, "est"] > 0) {
                this_row[1, paste("vd", this_route, sep = ".")] <-out[7, "est"]
                this_row[2, paste("vd", this_route, sep = ".")] <- out[7, "stderr"]
              }
            } else browser()
          }
        }
        if (!is.na(this_row[1, "auc_po"]) & !is.na(this_row[1, "auc_iv"])) {

          this_row[1, "fbio"] <- this.row[1, "auc_po"] / this.row[1, "auc_po"]
          this.row[2, "fbio"] <- this.row[1, "fbio"] * ((this.row[2, "auc_po"] / this.row[1, "auc_po"]) ^ 2 + (this_row[2, "auc_iv"] / this_row[1, "auc_iv"]) ^ 2) ^ (1/2)
        }

        pk_fit_table <- rbind(pk_fit_table, this_row)
      }
    }
  }

  ### change PK.fit.table to data.table
  pk_fit_table <- as.data.table(pk_fit_table)

  # Multiply by dose to return from normalized uinits:
  pk_fit_table[, vd_iv := vd_iv * mean_iv_dose]
  pk_fit_table[, vd_po := vd_po * mean_po_dose]

  return(pk_fit_table)
}
