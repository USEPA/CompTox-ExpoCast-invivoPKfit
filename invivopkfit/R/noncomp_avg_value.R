#' Helper function to take average value of replicate timepoints
#'
#' @param data A table of concentration-time data.
#'
#' @return A data.table with average concentration values
#'
#' @author Christopher Cook
#' @importFrom magrittr "%>%"

### function to take average of multiple timepoints
noncomp_avg_value <- function(data) {

  ### convert to data.frame, but fix this later because converting back and forth is silly
  data <- as.data.frame(data)

  ### if there are replicate timepoints, take the average Value for each timepoint
  if(length(unique(data$Time)) < length(data$Time)) {
    data <- data %>%
      dplyr::group_by(Time) %>%
      dplyr::do(dplyr::mutate(., mean = mean(.$Value, na.rm = TRUE)))

    ### if all values are NA, output will return 'NaN's
    ### coerces 'NaN's back to NA
    data$mean[is.nan(data$mean)] <- NA

    ### remove replicate timepoints
    data <- data %>%
      dplyr::distinct(Time, .keep_all = TRUE)

    ### delete Value column and reassign new mean column as Value column
    data$Value <- NULL
    names(data)[names(data) == "mean"] <- "Value"
  }

  data <- as.data.table(data)

  data
}

