estimate_conc_sd <- function(dat,
                         reference_col = "Reference",
                         chem_col = "DTXSID",
                         media_col = "Media",
                         species_col = "Species",
                         value_col = "Value",
                         sd_col = "Value_SD",
                         n_subj_col = "N_Subjects")
{

  #Impute missing SDs as the minimum non-missing SD per reference, chemical, species, and medium
  #the idea being that a missing SD probably occurred because there were no visible error bars on a point on a graph


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

  #for each group, impute any missing SD as the minimum non-missing SD.
  #do this only for observatiosn where N_Subjects > 1.
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
