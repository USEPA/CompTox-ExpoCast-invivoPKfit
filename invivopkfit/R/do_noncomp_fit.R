#' Helper function to perform noncomp fits
#'
#' @param data.set A table of concentration-time data.
#'
#' @return A data.table with average concentration values
#'
#' @author John Wambaugh, Caroline Ring, Christopher Cook
#' @importFrom magrittr "%>%"

do_noncomp_fit <- function(data.set) {

  ### for PK package to work, have to take average of replicate timepoints
  data.list <- list()
  for(this.compound in unique(data.set[, Compound])) {
    this.compound.data <- data.set[Compound == this.compound]

    ### list of columns by which to split data.set
    col_list <- list(this.compound.data$Reference,
                     this.compound.data$Dose,
                     this.compound.data$Route,
                     this.compound.data$Media)

    ### split data.set into list of data.tables
    split_list <- split(this.compound.data, col_list)

    ### remove data.tables with 0 rows
    # this.medium.list <- keep(this.medium.list, ~nrow(.) > 0)

    ### take average of replicate timepoints for each data.table
    split_list_avg <- lapply(split_list, noncomp_avg_value)
    sub_data <- do.call(rbind, split_list_avg)

    data.list[[this.compound]] <- sub_data
  }

  ### combine data.tables back into one large data.set
  data.set <- do.call(rbind, data.list)

  # Normalize all concentrations by dose to fit jointly:
  data.set[, Value.Norm := Value/Dose]

  ### not sure why PK.fit.table needs to be defined as NULL right here
  PK.fit.table <- NULL

  for (this.cas in sort(unique(data.set$CAS)))
  {
    this.subset <- subset(data.set, CAS == this.cas & !is.na(Value.Norm))
    for (this.species in sort(unique(this.subset$Species)))
    {
      # Need oral and iv to estimate bioavailability:
      if ("po" %in% unique(this.subset$Route) & "iv" %in% unique(this.subset$Route))
      {
        this.iv.time.list <- list()
        this.iv.conc.list <- list()
        this.iv.id.list <- list()
        this.oral.time.list <- list()
        this.oral.conc.list <- list()
        this.oral.id.list <- list()
        this.row <- data.frame(Compound = this.subset$Compound[1],
                               CAS = this.cas,
                               Reference = NA,
                               Source = NA,
                               Species = this.species,
                               Mean.iv.Dose = mean(subset(this.subset, Route == "iv")$Dose), ### averaged doses for iv and po
                               Mean.po.Dose = mean(subset(this.subset, Route == "po")$Dose),
                               model.type = "noncompartmental",
                               param.value.type = "Estimated",
                               AUC.po = NA,
                               Vd.po = NA,
                               CL.po = NA,
                               AUC.iv = NA,
                               Vd.iv = NA,
                               CL.iv = NA,
                               Fbio = NA,
                               stringsAsFactors = F)

        this.row <- rbind(this.row,this.row)

        rownames(this.row) <- NULL
        this.row[2,"param.value.type"] <- "Estimated std dev"
        for (this.route in unique(this.subset$Route))
        {
          this.route.subset <- subset(this.subset,Route == this.route)

          # new.route.subset <- NULL
          # print(this.route.subset$Compound)
          this.route.subset$Subject <- as.character(this.route.subset$Subject) ### added because Subject was logical and characters were not coerced correctly
          if (any(is.na(this.route.subset$Subject)))
          {
            for (this.reference in unique(this.route.subset$Reference))
            {
              for (this.dose in unique(this.route.subset$Dose))
              {
                this.route.subset[is.na(this.route.subset$Subject) & this.route.subset$Dose == this.dose, "Subject"] <- paste(this.reference, this.route, this.dose, sep = "-")
              }
            }
          }

          this.route.subset$id <- paste(this.row[1, "Reference"], this.route.subset$Subject, sep = "-")

          ### conc and time are required fields for the PK package
          this.route.subset$conc <- this.route.subset$Value.Norm
          this.route.subset$time <- this.route.subset$Time

          for (this.id in unique(this.route.subset$id))
          {
            ### comment here
            if (dim(subset(this.route.subset, id == this.id))[1] < 3)
            {
              this.route.subset <- subset(this.route.subset, id != this.id)
            }
          }

          ### order this.route.subset by time
          this.route.subset <- this.route.subset[order(this.route.subset$time), ]

          if (nrow(this.route.subset) > 0) {
            this.id.list <- list()

            for (this.id in unique(this.route.subset[, id]))
            {
              this.id.data <- this.route.subset[id == this.id]

              this.id.data <- as.data.frame(this.id.data)
              if (length(unique(this.id.data$time)) > 1) {
                this.id.data.group <- this.id.data %>%
                  dplyr::group_by(time) %>%
                  dplyr::tally()
                this.id.data <- dplyr::full_join(this.id.data, this.id.data.group)
                this.id.data <- as.data.table(this.id.data)
                this.id.data <- this.id.data[n == max(n)]
                this.id.data <- this.id.data[, n := NULL]
              }

              this.id.list[[this.id]] <- this.id.data
            }

            this.route.subset <- do.call(rbind, this.id.list)
          }

          if (dim(this.route.subset)[1] > 2)
          {
            if (this.route == "po")
            {
              dose.arg <- 0
            }
            else if (this.route == "iv")
            {
              dose.arg <- 1
            }

            ### depending on length time, try .complete or .batch
            if (length(this.route.subset$time) == length(unique(this.route.subset$time)))
            {
              out <- try(PK::nca.complete(data=this.route.subset,dose=dose.arg,method = "z"))
            }
            else
            {
              out <- try(PK::nca.batch(data = this.route.subset, dose = dose.arg, method = "z"))
            }


            if (class(out) != "try-error")
            {
              out <- out$CIs
              if (!is.na(out[1, "est"] > 0)) if (out[1, "est"] > 0)
              {
                this.row[1, paste("AUC", this.route, sep = ".")] <- out[1, "est"]
                this.row[2, paste("AUC", this.route, sep = ".")] <- out[1, "stderr"]
              }
              if (!is.na(out[6, "est"] > 0)) if (out[6, "est"] > 0)
              {
                this.row[1, paste("CL", this.route, sep = ".")] <- out[6, "est"]
                this.row[2, paste("CL", this.route, sep = ".")] <- out[6, "stderr"]
              }
              if(!is.na(out[7, "est"] > 0)) if (out[7, "est"] > 0)
              {
                this.row[1, paste("Vd", this.route, sep = ".")] <-out[7, "est"]
                this.row[2, paste("Vd", this.route, sep = ".")] <- out[7, "stderr"]
              }
            } else browser()
          }
        }
        if (!is.na(this.row[1, "AUC.po"]) & !is.na(this.row[1, "AUC.iv"]))
        {
          this.row[1, "Fbio"] <- this.row[1, "AUC.po"] / this.row[1, "AUC.iv"]
          this.row[2, "Fbio"] <- this.row[1, "Fbio"] * ((this.row[2, "AUC.po"] / this.row[1, "AUC.po"]) ^ 2 + (this.row[2, "AUC.iv"] / this.row[1, "AUC.iv"]) ^ 2) ^ (1/2)
        }

        ### if PK.fit.table is null until here, why define it earlier?
        PK.fit.table <- rbind(PK.fit.table,this.row)
      }
    }
  }

  ### change PK.fit.table to data.table
  PK.fit.table <- as.data.table(PK.fit.table)

  # Multiply by dose to return from normalized uinits:
  PK.fit.table[, Vd.iv := Vd.iv*Mean.iv.Dose]
  PK.fit.table[, Vd.po := Vd.po*Mean.po.Dose]

  return(PK.fit.table)
}
