#' Assign an estimated limit of quantification wherever missing
#'
#' Sometimes in mass spectrometry a compound is detected, but the signal is too
#' low and falls outside the range that can be interpreted quantitatively (that
#' is, as a concentration). These detects can still be used to refine parameter
#' estimates by preferring curve fits that predict concentrations below the
#' known "limit of quantitation" (LOQ) -- the lowest detectable quantifiable
#' concentration. LOQ's are not always available from data reported in the
#' literature, so this function estimates a LOQ such that all the reported
#' values are above the LOQ.

#' The LOQ is assumed to vary from reference to reference, as each study would
#' have different mass spectrometry equipment and methods. LOQ also is assumed
#' to vary by medium (e.g. blood or plasma), as preparation of samples may be
#' different. For each unique combination of reference, chemical, species, and
#' medium, the LOQ is estimated by multiplying the lowest detected concentration
#' (the minimum non-NA reported concentration greater than zero) by a constant
#' factor (\code{loq.factor}). If there are no detected concentrations for a
#' given combination of reference, chemical, species, and medium (i.e., all
#' reported concentrations are NA), then the imputed LOQ will be NA.
#'
#' @param dat A \code{data.frame} of concentration vs. time data with at least
#'   the following columns:
#'
#' @param reference_col The name of the column indicating different research
#'   study reference documents (and therefore different mass spectrometry
#'   methods). Default "document_reference.id".
#'
#' @param chem_col The name of the column indicating chemical identity. Default
#'   "chemicals_dosed.dsstox_substance_id".
#'
#' @param media_col The name of the column indicating sample medium. Default
#'   "series.conc_medium_normalized".
#'
#' @param value_col The name of the column indicating concentration. Default
#'   "conc_time_values.conc". Nondetects should be reported as NA in this
#'   column.
#'
#' @param species_col The name of the column indicating species. Default
#'   "subjects.species."
#'
#' @param loq_col The name of the LOQ column. Default "calc_loq". If this
#'   column does not already exist, it will be created.
#'
#' @param loq.factor A factor by which to multiply the lowest detected value in
#'   each group to estimate an LOQ. Default 0.45.
#'
#' @return A \code{data.frame} that is the same as \code{dat}, but with the
#'   column specifed in \code{calc_loq} added
#'
#' @author John Wambaugh
#'
#' @export estimate_loq
estimate_loq <- function(dat,
                         reference_col = "document_reference.id",
                         chem_col = "chemicals_dosed.dsstox_substance_id",
                         media_col = "series.conc_medium_normalized",
                         species_col = "subjects.species",
                         value_col = "conc_time_values.conc",
                         loq_col = "calc_loq",
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
  #detected concentration.
  dat_split_loq <- lapply(dat_split,
                     function(this_dat){
                       impute_loq <- calc_loq_factor *
                         min(this_dat[[value_col]], na.rm = TRUE)
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


  return(dat_out)
}
