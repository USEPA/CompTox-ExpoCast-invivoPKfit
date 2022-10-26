
rename_columns <- function(data.set,
                           compound.col,
                           cas.col,
                           reference.col,

                           species.col,
                           species.default,

                           species.weight.col,
                           species.weight.units.col,
                           species.weight.units.default,

                           dose.col,

                           time.col,
                           time.units.col,
                           time.units.default,

                           media.col,
                           media.units.col,
                           media.units.default,

                           value.col,

                           units.col,
                           units.default,

                           route.col,
                           route.default,

                           source.col,
                           source.default,

                           loq.col,
                           loq.default,

                           subject.col,
                           subject.default,

                           info.col,
                           info.default) {

  ### if default argument filled, create column and use default as value
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

  ### create vector of column names
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

  ### delete extra/NA columns
  data.set <- data.set[!is.na(names(data.set))]

  return(data.set)
}
