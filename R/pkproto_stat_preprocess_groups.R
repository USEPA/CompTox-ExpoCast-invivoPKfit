#' LOQ Group
#'
#' Defines the grouping variables for LOQ imputation in [do_preprocess.pk()].
#'
#' @param ...  A set of unquoted variables whose unique combinations define
#'  a group with which to impute missing LOQ values (using the
#'  minimum non-missing LOQ value in the group multiplied by `calc_loq_factor`).
#'  Default is `Chemical, Species, Reference, Media`.
#'
#' @returns A list of expressions. This is added to the `pk` object.
#' @export
#' @author Gilberto Padilla Mercado
#'
stat_loq_group <- function(...) {

  # get arguments and values
  this_group <- try(rlang::ensyms(...))
  if (inherits(this_group, "try-error") || rlang::is_empty(this_group)) {
    this_group <- rlang::syms(c("Chemical", "Species", "Reference", "Media"))
  }

  # set class
  class(this_group) <- c(class(this_group), "pkproto", "pk_loq_group")

  return(this_group)
}


#' SD Group
#'
#' Defines the grouping to calculate standard deviation of data in [do_preprocess.pk()].
#'
#' @param ... A set of unquoted variables whose unique combinations define
#'  a group with which to impute missing SD values (using the
#'  minimum non-missing SD value in the group).
#'  Default is `Chemical, Species, Reference, Media`.
#'
#' @returns A list of expressions. This is added to the `pk` object.
#' @export
#' @author Gilberto Padilla Mercado
#'
stat_sd_group <- function(...) {

  # get arguments and values
  this_group <- try(rlang::ensyms(...))
  if (inherits(this_group, "try-error") || rlang::is_empty(this_group)) {
    this_group <- rlang::syms(c("Chemical", "Species", "Reference", "Media"))
  }

  # set class
  class(this_group) <- c(class(this_group), "pkproto", "pk_sd_group")

  return(this_group)
}






