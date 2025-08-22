#' Get concentration scalings
#'
#' A helper function to get concentration scalings
#'
#' In methods applied to fitted [pk()] objects that also accept `newdata`
#' arguments, the user may specify whether to use the concentration scaling of
#' the fitted [pk()] object, or use a different concentration scaling. This is
#' done by specifying an argument `use_scale_conc`, which may be `TRUE` (to use
#' the scaling from the fitted object), `FALSE` (to use no scaling), or may be a
#' named list with elements `dose_norm` and `log10_trans` to specify
#' scaling/transformation directly. This helper function parses the
#' `use_scale_conc` argument.
#'
#' @param use_scale_conc The `use_scale_conc` argument (see Details)
#' @param obj A [pk()] object
#' @return A named list with elements `dose_norm` and `log10_trans`, both
#'   logical.
conc_scale_use <- function(use_scale_conc,
                           obj) {
  # apply scaling or not?
  if (isTRUE(use_scale_conc)) {
    dose_norm <- obj$scales$conc$dose_norm
    log10_trans <- obj$scales$conc$log10_trans
  } else if (isFALSE(use_scale_conc)) {
    dose_norm <- FALSE
    log10_trans <- FALSE
  } else {
    if (!is.list(use_scale_conc)) {
      cli::cli_abort(paste0(
        "{.var use_scale_conc} must be either TRUE, FALSE, or a list like ",
        "{.var list(dose_norm = TRUE/FALSE, log10_trans = TRUE/FALSE)`"
      ))
    }

    if (!all(c("dose_norm", "log10_trans") %in% names(use_scale_conc))) {
      cli_abort(paste0(
        "When {.var use_scale_conc} is a list, it must be named include both ",
           "{.val dose_norm} and {.val log10_trans} set to TRUE or FALSE."
        ))
    }

    if (!all(sapply(use_scale_conc[c("dose_norm", "log10_trans")], is.logical))) {
      cli_abort(paste0(
        "When {.var use_scale_conc} is a list with both {.val dose_norm} and {.val log10_trans}",
        "named elements, these must be set to TRUE or FALSE."
      ))
    }

    dose_norm <- use_scale_conc$dose_norm
    log10_trans <- use_scale_conc$log10_trans
  }

  return(list("dose_norm" = dose_norm,
              "log10_trans" = log10_trans))
}
