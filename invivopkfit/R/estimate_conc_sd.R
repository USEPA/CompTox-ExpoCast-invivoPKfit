#' Impute missing sample standard deviations
#'
#' For multi-subject summary data reported as sample mean, sample standard
#' deviation (SD), and number N of subjects, impute values for any missing
#' sample SDs.
#'
#' This function imputes missing SDs as the minimum non-missing SD per
#' reference, chemical, species, and medium. The idea is that a missing SD
#' probably occurred because there were no visible error bars on a point on a
#' graph, so data extractors/curators reported the SD missing. In this case, the
#' effective errorbars would be the size of the symbol itself. The minimum
#' non-missing SD is a conservative estimate of the symbol size; if error bars
#' were just barely visible, this would be the minimum non-missing SD.
#'
#' If there are no non-missing SDs for a given reference, chemical, species, and
#' medium, then this function imputes missing sample SDs as equal to the sample
#' means. This assumption was made after observation of a plot of sample SD vs.
#' sample mean for multi-subject data across chemicals and studies. The relation
#' is effectively on the identity line.
#'
#' If there are no non-missing SDs for a given reference, chemical, species, and
#' medium, and a sample mean is NA (missing or below LOD), then a corresponding
#' missing sample SD will remain missing.
#'
#' @param dat A `data.frame` containing the concentration vs. time data to be
#'   imputed. For example, from CvTdb.
#' @param reference_col Name of the variable in `dat` identifying unique
#'   references. Default "Reference".
#' @param chem_col Name of the variable in `dat` identifying unique
#'   chemicals/substances. Default "DTXSID."
#' @param media_col Name of the variable in `dat` identifying concentration
#'   media (blood or plasma). Default "Media".
#' @param species_col Name of the variable in `dat` identifying species. Default
#'   "Species".
#' @param value_col Name of the variable in `dat` that contains concentration
#'   values (sample means). Default "Value".
#' @param sd_col Name of the variable in `dat` that contains concentration
#'   sample standard deviation values. If this variable does not already exist
#'   in `dat`, it will be created and filled with `NA_real_`. Default
#'   "Value_SD".
#' @param n_subj_col Name of the variable in `dat` that contains the number of
#'   subjects for each sample mean/SD. Default "N_Subjects".
#' @return The same `data.frame` as `dat`, but with missing values in the
#'   variable `sd_col` imputed as described in Details.
#' @author Caroline Ring

estimate_conc_sd <- function(dat,
                         reference_col = "Reference",
                         chem_col = "DTXSID",
                         media_col = "Media",
                         species_col = "Species",
                         value_col = "Value",
                         sd_col = "Value_SD",
                         n_subj_col = "N_Subjects")
{




  #if SD column does not already exist, create it and fill with NA's
  if(!(sd_col %in% names(dat))){
    dat[[sd_col]] <- NA_real_
  }

  #split dat into groups by reference, chem, media, species
  dat_split <- split(dat,
                     dat[c(reference_col,
                           chem_col,
                           media_col,
                           species_col)],
                     drop = TRUE
  )

  #Impute missing SDs as the minimum non-missing SD per reference, chemical,
  #species, and medium. The idea is that a missing SD probably occurred because
  #there were no visible error bars on a point on a graph. In this case, the
  #effective errorbars would be the size of the symbol itself. The minimum
  #non-missing SD is a conservative estimate of the symbol size; if error bars
  #were just barely visible, this would be the minimum non-missing SD.

  dat_split_sd <- lapply(dat_split,
                          function(this_dat){
                            suppressWarnings(impute_sd <- min(this_dat[[sd_col]], na.rm = TRUE))
                            #If there were no detects for this group, the imputed SD
                            #will be Inf. Set it to NA instead.
                            if(!is.finite(impute_sd)){
                              impute_sd <- NA_real_
                            }
                            this_dat[(this_dat[[n_subj_col]] > 1) %in% TRUE &
                              is.na(this_dat[[sd_col]]),
                                     sd_col] <- impute_sd
                            return(this_dat)
                          })

  #rejoin
  dat_out <- unsplit(dat_split_sd,
                     dat[c(reference_col,
                           chem_col,
                           media_col,
                           species_col)],
                     drop = TRUE)

  #For any remaining cases with N_Subjects > 1 and SD still NA, impute SD as
  #mean itself. This is based on examining the relationship of mean and SD for
  #individual animal data by reference, chemical, dose, and time. The
  #relationship is nearly the identity line.
if(any((dat_out[[n_subj_col]] > 1) %in% TRUE &
       is.na(dat_out[[sd_col]]))){
  selected <- (dat_out[[n_subj_col]] > 1) %in% TRUE &
    is.na(dat_out[[sd_col]])
  dat_out[selected, sd_col] <- dat_out[selected, value_col]
}

  return(dat_out)
}
