#' Helper function to calculate LOQ for processed data
#'
#' @param data data.table with Value and LOQ column names
#'
#' @return A data.table with assigned LOQ values for those that were NA
#'
#' @author Christopher Cook

fix_loq <- function(data) {

  data <- as.data.frame(data)

  if (all(is.na(data$LOQ))) {data$LOQ <- 0.45 * min(data$Value, na.rm = TRUE)}

  data <- as.data.table(data)

  return(data)
}
