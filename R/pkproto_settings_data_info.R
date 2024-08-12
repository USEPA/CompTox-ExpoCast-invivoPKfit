#' Data info settings
#'
#' @param nca_group A set of variables. Data will be split into groups according
#'   to unique combinations of these variables, and non-compartmental analysis
#'   will be performed separately on each group. Default `dplyr::vars(Chemical,
#'   Species, Reference, Route, Media, Dose)`.
#' @param ... Other arguments (currently ignored).
#' @return An object of class `c(pkproto, pk_settings_data_info)`
#' @export
settings_data_info <- function(nca_group = dplyr::vars(Chemical,
                                                       Species,
                                                       Reference,
                                                       Route,
                                                       Media,
                                                       Dose),
                               summary_group = dplyr::vars(Chemical,
                                                           Species,
                                                           Route,
                                                           Media,
                                                           Dose),
                               ...){
  Chemical <- NULL
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_settings_data_info <- argg
  #set class
  class(this_settings_data_info) <- c(class(this_settings_data_info),
                                      "pkproto",
                                      "pk_settings_data_info")

  return(this_settings_data_info)
}

