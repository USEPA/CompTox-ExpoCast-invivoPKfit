
rename_columns <- function(data,

                           compound_col = compound_col,
                           compound_default = compound_default,

                           cas_col = cas_col,
                           cas_default = cas_default,

                           dtxsid_col = dtxsid_col,
                           dtxsid_default = dtxsid_default,

                           reference_col = reference_col,
                           reference_default = reference_default,

                           species_col = species_col,
                           species_default = species_default,

                           species_weight_col = species_weight_col,
                           species_weight_default = species_weight_default,
                           species_weight_units_col = species_weight_units_col,
                           species_weight_units_default = species_weight_units_default,

                           dose_col = dose_col,
                           dose_default = dose_default,

                           time_col = time_col,
                           time_default = time_default,
                           time_units_col = time_units_col,
                           time_units_default = time_units_default,

                           media_col = media_col,
                           media_default = media_default,
                           media_units_col = media_units_col,
                           media_units_default = media_units_default,

                           conc_col = conc_col,
                           conc_default = conc_default,
                           conc_units_col = conc_units_col,
                           conc_units_default = conc_units_default,

                           route_col = route_col,
                           route_default = route_default,

                           loq_col = loq_col,
                           loq_default = loq_default,

                           loq_units_col = loq_units_col,
                           loq_units_default = loq_units_default) {

  ### if default argument filled, create column and use default as value

  ### compound
  if(!is.null(compound_default)) {
    data[, compound_col] <- compound_default
  }

  ### cas
  if(!is.null(cas_default)) {
    data[, cas_col] <- cas_default
  }

  ### reference
  if(!is.null(reference_default)) {
    data[, reference_col] <- reference_default
  }

  ### species
  if(!is.null(species_default)) {
    data[, species_col] <- species_default
  }

  ### species_weight_units
  if(!is.null(species_weight_units_default)) {
    data[, species_weight_units_col] <- species_weight_units_default
  }

  ### time
  if(!is.null(time_default)) {
    data[, time_col] <- time_default
  }

  ### time_units
  if(!is.null(time_units_default)) {
    data[, time_units_col] <- time_units_default
  }

  ### media
  if(!is.null(media_units_default)) {
    data[, media_units_col] <- media_units_default
  }

  ### media
  if(!is.null(media_default)) {
    data[, media_col] <- media_default
  }

  ### conc_units
  if(!is.null(conc_default)) {
    data[, conc_col] <- conc_default
  }

  ### conc_units
  if(!is.null(conc_units_default)) {
    data[, conc_units_col] <- conc_units_default
  }

  ### route
  if(!is.null(route_default)) {
    data[, route_col] <- route_default
  }

  ### loq
  if(!is.null(loq_default)) {
    data[, loq_col] <- loq_default
  }

  ### loq_units
  if(!is.null(loq_units_default)) {
    data[, loq_units_col] <- loq_units_default
  }

  ### create vector of column names
  cols <- c(compound_col,
            cas_col,
            reference_col,
            species_col,
            species_weight_col,
            species_weight_units_col,
            dose_col,
            time_col,
            time_units_col,
            media_col,
            media_units_col,
            conc_col,
            conc_units_col,
            route_col,
            loq_col,
            loq_units_col)

  ### stop if any columns missing from data.set
  if (!(all(cols %in% colnames(data)))) {
    stop(paste("Missing columns named:",
               paste(cols[!(cols %in% colnames(data))], collapse = ", ")))
  }

  ### Set column order of data.table to cols vector
  data <- data.table::setcolorder(data, cols)

  # Standardize the column names:
  compound_col <- "compound"
  cas_col <- "cas"
  reference_col <- "reference"
  species_col <- "species"
  species_weight_col <- "species_weight"
  species_weight_units_col <- "species_weight_units"
  dose_col <- "dose"
  time_col <- "time"
  time_units_col <- "time_units"
  media_col <- "media"
  media_units_col <- "media_units"
  conc_col <- "conc"
  conc_units_col <- "conc_units"
  route_col <- "route"
  loq_col <- "loq"
  loq_units_col <- "loq_units"

  ### rename colnames to standardized colnames
  colnames(data) <- c(compound_col,
                      cas_col,
                      reference_col,
                      species_col,
                      species_weight_col,
                      species_weight_units_col,
                      dose_col,
                      time_col,
                      time_units_col,
                      media_col,
                      media_units_col,
                      conc_col,
                      conc_units_col,
                      route_col,
                      loq_col,
                      loq_units_col)

  ### delete extra/NA columns
  data <- data[!is.na(names(data))]

  return(data)
}
