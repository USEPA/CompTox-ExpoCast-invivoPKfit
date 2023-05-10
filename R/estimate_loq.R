#' Assign an estimated limit of quantification wherever missing
#'
#' Sometimes in mass spectrometry a compound is detected, but the signal is too
#' low and falls outside the range that can be interpreted quantitatively (that
#' is, as a concentration). These detects can still be used to refine parameter
#' estimates by preferring curve fits that predict concentrations below the
#' known "limit of quantitation" (LOQ) -- the lowest detectable quantifiable
#' concentration. LOQs are not always available from data reported in the
#' literature. When LOQ is missing, this function imputes an LOQ such that all
#' the reported values are above the LOQ.

#' The LOQ is assumed to vary from reference to reference, as each study would
#' have different mass spectrometry equipment and methods. LOQ also is assumed
#' to vary by medium (e.g. blood or plasma), as preparation of samples may be
#' different for different media.
#'
#' For each unique combination of reference, chemical, species, and medium, the
#' LOQ is estimated by multiplying the lowest detected concentration (the
#' minimum non-NA reported concentration greater than zero) by a constant factor
#' (`calc_loq_factor`). If there are no detected concentrations for a given
#' combination of reference, chemical, species, and medium (*i.e.*, all reported
#' concentrations are NA), then the imputed LOQ will be NA.
#'
#' @param dat A \code{data.frame} of concentration vs. time data with at least
#'   the following variables:
#'
#' @param reference_col The name of the variable in `dat` indicating different
#'   research study reference documents (and therefore different mass
#'   spectrometry methods). Default "Reference".
#'
#' @param chem_col The name of the variable in `dat` indicating chemical
#'   identity. Default "DTXSID".
#'
#' @param media_col The name of the variable in `dat` indicating sample medium.
#'   Default "Medium".
#'
#' @param value_col The name of the variable in `dat` indicating concentration.
#'   Default "Value". Nondetects should be reported as NA in this column.
#'
#' @param species_col The name of the variable in `dat` indicating species.
#'   Default "Species".
#'
#' @param loq_col The name of the LOQ column. Default "LOQ". If this column does
#'   not already exist, it will be created and filled with `NA_real_`.
#'
#' @param calc_loq_factor A factor by which to multiply the lowest detected
#'   value in each group to estimate an LOQ. Default 0.45.
#'
#' @return A \code{data.frame} that is the same as \code{dat}, but with any
#'   missing values in variable `calc_loq` imputed as above.
#'
#' @author John Wambaugh, Caroline Ring
#'
#' @export estimate_loq
estimate_loq <- function(dat,
                         reference_col = "Reference",
                         chem_col = "DTXSID",
                         media_col = "Media",
                         species_col = "Species",
                         value_col = "Value",
                         loq_col = "LOQ",
                         calc_loq_factor = 0.45)
{

  #if LOQ column does not already exist, create it and fill with NA's
  if(!(loq_col %in% names(dat))){
    dat[[loq_col]] <- NA_real_
  }

  #split dat into groups by reference, chem, media, species
  dat_split <- split(dat,
                     dat[c(reference_col,
                           chem_col,
                           media_col,
                           species_col)],
                     drop = TRUE
                     )

  #for each group, impute any missing LOQ as calc_loq_factor * the minimum
  #detected concentration (that is greater than zero).
  dat_split_loq <- lapply(dat_split,
                     function(this_dat){
                       these_values <- this_dat[[value_col]]
                       suppressWarnings(
                         minval <- min(these_values[these_values > 0],
                                       na.rm = TRUE))
                       impute_loq <- calc_loq_factor * minval
                       #If there were no detects for this group, the imputed LOQ
                       #will be Inf. Set it to NA instead.
                       if(!is.finite(impute_loq)){
                         impute_loq <- NA_real_
                       }
                       this_dat[is.na(this_dat[[loq_col]]),
                                loq_col] <- impute_loq
                       return(this_dat)
                     })

  #rejoin
  dat_out <- unsplit(dat_split_loq,
                     dat[c(reference_col,
                           chem_col,
                           media_col,
                           species_col)],
                     drop = TRUE)

  #If all observations were NA for a given reference, chemical, media, and
  #species, then imputed LOQ will also be NA.

  return(dat_out)
}
