get_starts_flat <- function(data,
                            par_DF){

  Vdist <- NA_real_
  Fgutabs_Vdist <- NA_real_
  Fgutabs <- NA_real_
  Rblood2plasma <- 1

  #Work only with detects for these rough estimates
  tmpdat <- subset(data,
                   Detect %in% TRUE)

  #Split into IV and PO
  ivdat <- subset(tmpdat,
                  Route %in% "iv")
  podat <- subset(tmpdat,
                  Route %in% "oral")

  has_iv <- any(tmpdat$Route %in% "iv")
  has_po <- any(tmpdat$Route %in% "oral")
  has_plasma <- any(tmpdat$Media %in% "plasma")
  has_blood <- any(tmpdat$Media %in% "blood")
  has_iv_plasma <- any(tmpdat$Media %in% "plasma" &
                      tmpdat$Route %in% "iv")
  has_iv_blood <- any(tmpdat$Media %in% "blood" &
                         tmpdat$Route %in% "iv")
  has_po_plasma <- any(tmpdat$Media %in% "plasma" &
                         tmpdat$Route %in% "oral")
  has_po_blood <- any(tmpdat$Media %in% "blood" &
                        tmpdat$Route %in% "oral")

  Vdist_plasma_log10 <- NA_real_
  Vdist_blood_log10 <- NA_real_

  Fgutabs_Vdist_plasma_log10 <- NA_real_
  Fgutabs_Vdist_blood_log10 <- NA_real_

  Rblood2plasma_iv_log10 <- NA_real_
  Rblood2plasma_po_log10 <- NA_real_

  if(has_iv %in% TRUE){
    if(has_iv_plasma %in% TRUE){
      #get Vdist for plasma
      Vdist_plasma_log10 <- with(subset(ivdat,
                               Media %in% "plasma"),
                        -1 * mean(log10(Conc/Dose)),
                             na.rm = TRUE)
    }
      if(has_iv_blood %in% TRUE){ #if both blood and plasma IV data, estimate Rblood2plasma
        #get Vdist for blood
        Vdist_blood_log10 <- with(subset(ivdat,
                             Media %in% "blood"),
                      -1 * mean(log10(Conc/Dose)),
                      na.rm = TRUE)
      }
    if(has_iv_plasma %in% TRUE &
       has_iv_blood %in% TRUE){
      Rblood2plasma_iv_log10 <- Vdist_plasma_log10 + Vdist_blood_log10
    }
  }

  if(has_po %in% TRUE){
    if(has_po_plasma %in% TRUE){
      #get Fgutabs_Vdist for plasma
      Fgutabs_Vdist_plasma_log10 <- with(subset(podat,
                                        Media %in% "plasma"),
                                 mean(log10(Conc/Dose)),
                                 na.rm = TRUE)
    }
    if(has_po_blood %in% TRUE){
      #get Fgutabs_Vdist for blood
      Fgutabs_Vdist_blood_log10 <- with(subset(podat,
                                       Media %in% "blood"),
                                mean(log10(Conc/Dose)),
                                na.rm = TRUE)
    }
    if(has_po_plasma %in% TRUE &
       has_po_blood %in% TRUE){
      Rblood2plasma_po_log10 <- Fgutabs_Vdist_blood_log10 + Fgutabs_Vdist_blood_log10
    }
  }

if(has_iv %in% TRUE){
      Vdist <- 10^(mean(c(Vdist_plasma_log10,
                        Vdist_blood_log10),
                        na.rm = TRUE))
}

if(has_po %in% TRUE){
  Fgutabs_Vdist <- 10^(mean(c(Fgutabs_Vdist_plasma_log10,
                      Fgutabs_Vdist_blood_log10),
                      na.rm = TRUE))
}

  if(has_iv %in% TRUE &
     has_po %in% TRUE){
    Fgutabs <- Fgutabs_Vdist * Vdist
  }

  if(has_plasma %in% TRUE &
     has_blood %in% TRUE){
  Rblood2plasma <- 10^(mean(c(Rblood2plasma_iv_log10,
                              Rblood2plasma_po_log10),
                            na.rm = TRUE))
  }

  #update starting Concs
par_DF["Vdist", "start"] <- Vdist
par_DF["Fgutabs_Vdist", "start"] <- Fgutabs_Vdist
par_DF["Fgutabs", "start"] <- Fgutabs
par_DF["Rblood2plasma", "start"] <- Rblood2plasma

return(par_DF)
}
