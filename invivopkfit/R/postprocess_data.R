postprocess_data <- function(PK.fit.bind,
                             data.set){
#Do post fit calculations:
### basically calculate parameters based on previously calculated paramaters
PK.fit.table <- PK.fit.bind
if (model == "1compartment") {
  #Get stats for fitted total clearance :    L/kg body weight/h
  PK.fit.table[, CLtot := signif(Vdist * kelim, sig.figs)]

  #      PK.fit.table[,
  #                   AUC1mgkg:=kgutabs/Vdist/(kgutabs-kelim)*(
  #                     exp(-kgutabs*Tmax)/kgutabs  -exp(-kelim*Tmax)/kelim-
  #                     1/kgutabs + 1/kelim)]
  #      if (!is.na(Fgutabs)) PK.fit.table[,AUC1mgkg:=Fgutabs*AUC1mgkg]


  #Get statistics for Css from fitted CLtot values
  #Get stats for fitted total clearance:
  ### if Fgutabs is NA, assign value of 1, else return existing Fgutabs value, then divide by (24 * CLtot)
  PK.fit.table[,
               Css := ifelse(is.na(Fgutabs), 1, Fgutabs) / (24 * CLtot)] # 1 mg/kg/day / L/day/kg -> mg/L

  #Get statistics for halflife from fitted values
  PK.fit.table[,
               halflife := signif(log(2) / kelim, sig.figs)]


  PK.fit.table[,
               tpeak.oral := signif(log(kgutabs / kelim) / (kgutabs - kelim), sig.figs)]

  ### why does Fgutabs have to be 1 right here
  ### this just copies the Ccompartment column
  ### not sure how values correspond to specific param.value.types
  ### again, it's just a time, conc, and auc matrix
  PK.fit.table[,
               Cpeak.oral.1mgkg := signif(analytic_1comp_fun(
                 params=list(Fgutabs = 1,
                             kgutabs = kgutabs,
                             kelim = kelim,
                             Vdist = Vdist),
                 dose = 1,
                 tpeak.oral,
                 iv.dose = F)[, "Ccompartment"], sig.figs)]

  ### subset only certain param.value.types, not sure why
  PK.fit.table <- PK.fit.table[param.value.type %in% c("Predicted",
                                                       "Fitted geometric mean",
                                                       "Fitted geometric std dev")]

} else if(model == "2compartment") {
  PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                       "Fitted geometric mean",
                                       "Fitted mode"),
               c("beta",
                 "alpha") := lapply(list(Fbetaofalpha*Ralphatokelim*kelim,
                                         Ralphatokelim*(kelim+10^-6)),
                                    function(x) signif(x, sig.figs))]

  PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                       "Fitted geometric mean",
                                       "Fitted mode"),
               c("k21",
                 "k12") := lapply(list(alpha * beta / kelim,
                                       alpha + beta - kelim - alpha * beta / kelim),
                                  function(x) signif(x, sig.figs))]

  PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                       "Fitted geometric mean",
                                       "Fitted mode"),
               c("halflife",
                 "Vss",
                 "CLtot",
                 "Varea.or.Vbeta") := lapply(list(log(2) / beta,
                                                  V1 * (k21 + k12) / k21,
                                                  V1 * (k21 + k12) / k21 * beta,
                                                  V1 * kelim / beta),
                                             function(x) signif(x, sig.figs))]
  #Get Css = average plasma concentration for 1 mg/kg/day every 1 days
  #this is the same as the average for an equivalent constant oral infusion --
  #(1/24) mg/kg/hour every hour,
  #(1/(24*60)) mg/kg/minute every minute,
  #(1/(24*60*60)) mg/kg/second every second,
  #however finely you want to subdivide it.
  #If you subdivide it finely enough then you won't get peaks and valleys,
  #you'll just get an overall average time course.
  #the 1 in the numerator = dose = 1 mg/kg/day
  #the 1 in the denominator = time interval between doses = 1 day
  PK.fit.table[, Css := ifelse(is.na(Fgutabs), 1, Fgutabs) / (24 * CLtot)]

  PK.fit.table <- PK.fit.table[param.value.type %in% c("Predicted",
                                                       "Fitted geometric mean",
                                                       "Fitted geometric std dev")]
} else if (model == 'flat') {
  PK.fit.table <- PK.fit.table[param.value.type %in% c("Predicted",
                                                       "Fitted geometric mean",
                                                       "Fitted geometric std dev")]
}

###########################################################################################################
###########################################################################################################
###########################################################################################################

PK.fit.table <- PK.fit.table[order(Compound, Species)]

### split data.set into list of data.frame objects
split_df <- split(data.set, list(data.set$Compound, data.set$Reference, data.set$Media), drop = TRUE)

### apply fix_loq to each data.set
split_df_loq <- lapply(split_df, fix_loq)

### 'unsplit' data.sets
data.set <- do.call(rbind, split_df_loq)

rownames(data.set) <- c()

### coerce data.set back to data.table
# data.set <- as.data.table(data.set)

out <- list(PK.fit.table, data.set)

return(out)
}
