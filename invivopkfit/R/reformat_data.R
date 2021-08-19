reformat_data <- function(data.set,

                          compound.col = "Compound",
                          cas.col = "CAS",
                          reference.col = "Reference",

                          species.col = "Species",
                          species.default = NULL,

                          species.weight.col = "Species.Weight",
                          species.weight.units.col = "Species.Weight.Units",
                          species.weight.units.default = NULL,

                          dose.col = "Dose",

                          time.col = "Time",
                          time.units.col = "Time.Units",
                          time.units.default = NULL,

                          media.col = "Media",
                          media.units.col = "Media.Units",
                          media.units.default = NULL,

                          value.col = "Value",

                          units.col = "Units",
                          units.default = NULL,

                          route.col = "Route",
                          route.default = NULL,

                          source.col = "Source",
                          source.default = NULL,

                          loq.col = "LOQ",
                          loq.default = NULL,

                          subject.col = "Subject",
                          subject.default = NULL,

                          info.col = "info",
                          info.default = NULL) {

  if(!is.null(species.default)) {
    data.set[, species.col] <- species.default
  }

  if(!is.null(species.weight.units.default)) {
    data.set[, species.weight.units.col] <- species.weight.units.default
  }

  if(!is.null(time.units.default)) {
    data.set[, time.units.col] <- time.units.default
  }

  if(!is.null(media.units.default)) {
    data.set[, media.units.col] <- media.units.default
  }

  if(!is.null(units.default)) {
    data.set[, units.col] <- units.default
  }

  if(!is.null(route.default)) {
    data.set[, route.col] <- route.default
  }

  if(!is.null(loq.default)) {
    data.set[, loq.col] <- loq.default
  }

  if(!is.null(subject.default)) {
    data.set[, subject.col] <- subject.default
  }

  if(!is.null(info.default)) {
    data.set[, info.col] <- info.default
  }

  cols <- c(
    compound.col,
    cas.col,
    reference.col,
    species.col,
    species.weight.col,
    species.weight.units.col,
    dose.col,
    time.col,
    time.units.col,
    media.col,
    media.units.col,
    value.col,
    units.col,
    route.col,
    source.col,
    loq.col,
    subject.col,
    info.col
  )

  ### stop if any columns missing from data.set
  if (!(all(cols %in% colnames(data.set))))
  {
    stop(paste("Missing columns named:",
               paste(cols[!(cols%in%colnames(data.set))],collapse=", ")))
  }

  ### Set column order of data.table to cols vector
  data.set <- data.table::setcolorder(data.set, cols)

  # Standardize the column names:
  compound.col <- "Compound"
  cas.col <- "CAS"
  reference.col <- "Reference"
  species.col <- "Species"
  species.weight.col <- "Species.Weight"
  species.weight.units.col <- "Species.Weight.Units"
  dose.col <- "Dose"
  time.col <- "Time"
  time.units.col <- "Time.Units"
  media.col <- "Media"
  media.units.col <- "Media.Units"
  value.col <- "Value"
  units.col <- "Units"
  route.col <- "Route"
  source.col <- "Source"
  loq.col <- "LOQ"
  subject.col <- "Subject"
  info.col <- "info"

  ### rename colnames to standardized colnames
  colnames(data.set) <- c(
    compound.col,
    cas.col,
    reference.col,
    species.col,
    species.weight.col,
    species.weight.units.col,
    dose.col,
    time.col,
    time.units.col,
    media.col,
    media.units.col,
    value.col,
    units.col,
    route.col,
    source.col,
    loq.col,
    subject.col,
    info.col
  )

  ### rename extra/NA columns in template as "delete"
  names_corrected <- tidyr::replace_na(colnames(data.set), "delete")
  colnames(data.set) <- names_corrected

  ### delete columns named "delete"
  if("delete" %in% colnames(data.set)) {data.set[, "delete" == names(data.set)] <- NULL}

  ### Coerce all 'Value' values to be numeric
  data.set$Value <- as.numeric(data.set$Value)

  data.set <- as.data.table(data.set)

  return(data.set)
}
