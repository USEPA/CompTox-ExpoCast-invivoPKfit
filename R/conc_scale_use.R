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
      stop("use_scale_conc must be either TRUE, FALSE, ",
           "or a list like `list(dose_norm = TRUE/FALSE, ",
           "log10_trans = TRUE/FALSE)`")
    }

    if (!all(c("dose_norm", "log10_trans") %in% names(use_scale_conc))) {
      stop("use_scale_conc must be either TRUE, FALSE, ",
           "or a list like `list(dose_norm = TRUE/FALSE, ",
           "log10_trans = TRUE/FALSE)`")
    }

    if (!all(sapply(use_scale_conc[c("dose_norm", "log10_trans")], is.logical))) {
      stop("use_scale_conc must be either TRUE, FALSE, ",
           "or a list like `list(dose_norm = TRUE/FALSE, ",
           "log10_trans = TRUE/FALSE)`")
    }

    dose_norm <- use_scale_conc$dose_norm
    log10_trans <- use_scale_conc$log10_trans
  }

  return(list("dose_norm" = dose_norm,
              "log10_trans" = log10_trans))
}
