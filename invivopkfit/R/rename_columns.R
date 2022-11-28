#' Produce data with harmonized variable names.
#'
#' Renames data to use harmonized variable names.
#'
#' Variables will be renamed according to the following table:
#'
#' | Old name | New name | Definition |
#' | -------- | -------- | ----------------------|
#' | `compound.col` | "Compound" | Compound name |
#' | `cas.col` | "CAS" | CASRN |
#' | `reference.col` | "Reference" | Reference identifier |
#' | `species.col` | "Species" | Species name |
#' | `species.weight.col` | "Species.Weight" | Species-specific body weight |
#' | `species.weight.units.col` | "Species.Weight.Units" | Units of species-specific body weight |
#' | `dose.col` | "Dose" | Administered dose |
#' | `time.col` | "Time" | Time of each observation |
#' | `time.units.col` | "Time.Units" | Units of times |
#' |`media.col` | "Media" | Medium in which concentration was observed (e.g. "plasma" or "blood")
#' | `media.units.col` | "Media.Units" | Units of media (e.g. "harmonized") |
#' | `value.col` | "Value" | Observed concentrations |
#' | `units.col` | "Units" | Units of observed concentrations |
#' | `route.col` | "Route" | Route of dose administration (e.g. "iv" or "po") |
#' | `source.col` | "Source" | Source identifier (may be the same as reference identifier) |
#' | `loq.col` | "LOQ" | Limit of quantification for each observed concentration |
#' | `subject.col` | "Subject" | Subject identifier, if available (e.g. individual animal ID) |
#' | `info.col` | "Info" | Any other information for each observation |
#'
#' @param data.set A \code{data.frame} of concentration-time data.
#' @param compound.col Column name in \code{data.set} that identifies chemical
#'   compound. Default "Compound".
#' @param cas.col Column name in \code{data.set} to identify CASRN.
#'   Default"CAS".
#' @param reference.col Column name in \code{data.set} to identify reference.
#'   Default "Reference".
#' @param species.col Column name in \code{data.set} to identify species.
#'   Default "Species".
#' @param species.default If no species column exists in \code{data.set}, one
#'   will be created and filled with this value.  Default NULL.
#' @param species.weight.col Column name in \code{data.set} to identify species
#'   weight. Default "Species.Weight".
#' @param species.weight.units.col Column name in \code{data.set} to identify
#'   species weigh units. Default"Species.Weight.Units".
#' @param species.weight.units.default If no species weight units column exists
#'   in \code{data.set}, one will be created and filled with this value. Default
#'   NULL.
#' @param dose.col Column name in \code{data.set} to identify dose. Default
#'   "Dose".
#' @param time.col Column name in \code{data.set} to identify time. Default
#'   "Time".
#' @param time.units.col Column name in \code{data.set} to identify time units.
#'   Default "Time.Units."
#' @param time.units.default If no time units column exists in \code{data.set},
#'   one will be created and filled with this value.  Default NULL.
#' @param media.col Column name in \code{data.set} to identify media. Default
#'   "Media".
#' @param media.units.col Column name in \code{data.set} to identify media
#'   units. Default "Media.Units".
#' @param media.units.default If no media units column exists in
#'   \code{data.set}, one will be created and filled with this value.  Default
#'   NULL.
#' @param value.col Column name in \code{data.set} to identify value. Default
#'   "Value".
#' @param units.col Column name in \code{data.set} to identify units. Default
#'   "Units".
#' @param units.default If no units column exists in \code{data.set}, one will
#'   be created and filled with this value.  Default NULL.
#' @param route.col Column name in \code{data.set} to identify route of
#'   administration. Default "Route".
#' @param route.default If no route column exists in \code{data.set}, one will
#'   be created and filled with this value.  Default NULL.
#' @param source.col Column name in \code{data.set} to identify source. Default
#'   "Source."
#' @param source.default If no source column exists in \code{data.set}, one will
#'   be created and filled with this value.  Default NULL.
#' @param loq.col Column name in \code{data.set} to identify LOQ. Default "LOQ".
#' @param loq.default If no LOQ column exists in \code{data.set}, one will be
#'   created and filled with this value.  Default NULL.
#' @param subject.col Column name in \code{data.set} to identify subject.
#'   Default "Subject."
#' @param subject.default If no subject column exists in \code{data.set}, one
#'   will be created and filled with this value.  Default NULL.
#' @param info.col Column name in \code{data.set} to serve as info column.
#'   Default "Info".
#' @param info.default If no info column exists in \code{data.set}, one will be
#'   created and filled with this value.  Default NULL.
#'
#' @return A `data.frame` which is `data.set` with variables renamed to
#'   harmonized names.
#' @author John Wambaugh, Chris Cook

rename_columns <- function(data.set,
                           compound.col,
                           dtxsid.col,
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

  data.set <- as.data.frame(data.set)

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
    dtxsid.col,
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
  data.set <- data.set[, cols]
  #data.set <- data.table::setcolorder(data.set, cols)

  # Standardize the column names:
  compound.col <- "Compound"
  dtxsid.col <- "DTXSID"
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
    dtxsid.col,
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
