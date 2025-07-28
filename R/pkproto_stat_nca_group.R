#' NCA group
#'
#' @param ... A set of variables. Data will be split into groups according
#'   to unique combinations of these variables, and non-compartmental analysis
#'   will be performed separately on each group. Default `Chemical,
#'   Species, Reference, Route, Media, Dose`.
#' @return An object of class `c(pkproto, pk_settings_data_info)`
#' @export
stat_nca_group <- function(...) {
  # get arguments and values
  this_group <- try(rlang::ensyms(...))

  if (inherits(this_group, "try-error") || rlang::is_empty(this_group)) {
    this_group <- rlang::syms(c("Chemical", "Species", "Reference", "Route", "Media", "Dose"))
  }

  # set class
  class(this_group) <- c(class(this_group),
                         "pkproto",
                         "pk_nca_group")

  return(this_group)
}
