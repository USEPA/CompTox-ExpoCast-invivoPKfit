#' Main fitting function
#'
#' Fits parameters of a specified model to concentration-time data given in
#' data.set
#'
#' @param data.set A table of concentration-time data. Preferably
#'   \code{pkdataset_nheerlcleaned} or \code{pkdataset_nheerlorig}.
#' @param model Which general model should be fit for each chemical. Presently,
#'   only "1compartment" and "2compartment" are implemented.
#' @param modelfun Either "analytic" or "full" -- whether to fit using the
#'   analytic solution to the model, or the full ODE model. Presently,
#'   "analytic" is recommended (because the analytic solution is exact and much
#'   faster).
#' @param ratio.data.to.dose between the weights used to report the data and the weights used
#'   for the dose. For example, ug/L data and mg/kg/day dose would be 0.001
#'   (defaults to 1)
#'
#' @param compound.col Column in data.set to serve as compound column
#' @param cas.col Column in data.set to serve as CAS column
#' @param reference.col Column in data.set to serve as reference column
#' @param species.col Column in data.set to serve as species column
#' @param species.default Default value to fill newly created species column, should data.set lack a column equivalent to species
#' @param species.weight.col Column in data.set to serve as species.weight column
#' @param species.weight.units.col Column in data.set to serve as species.weight.units column
#' @param species.weight.units.default Default value to fill newly created species.weight.units column, should data.set lack a column equivalent to species.weight.units
#' @param dose.col Column in data.set to serve as dose column
#' @param time.col Column in data.set to serve as time column
#' @param time.units.col Column in data.set to serve as time.units column
#' @param time.units.default Default value to fill newly created time.units column, should data.set lack a column equivalent to time.units
#' @param media.col Column in data.set to serve as media column
#' @param media.units.col Column in data.set to serve as media.units column
#' @param media.units.default Default value to fill newly created media.units coloumn, should data.set lack a column equivalent to media.units
#' @param value.col Column in data.set to serve as value column
#' @param units.col Column in data.set to serve as units column
#' @param units.default Default value to fill newly created units column, should data.set lack a column equivalent to units
#' @param route.col Column in data.set to serve as route column
#' @param route.default Default value to fill newly created route column, should data.set lack a column equivalent to route
#' @param source.col Column in data.set to serve as source column
#' @param source.default Default value to fill newly created source column, should data.set lack a column equivalent to source
#' @param loq.col Column in data.set to serve as loq column
#' @param loq.default Default value to fill newly created loq column, should data.set lack a column equivalent to loq
#' @param subject.col Column in data.set to serve as subject column
#' @param subject.default Default value to fill newly created subject column, should data.set lack a column equivalent to subject
#' @param info.col Column in data.set to serve as info column
#' @param info.default Default value to fill newly created info column, should data.set lack a column equivalent to info
#'
#' @return A data.table of fitted parameter values for each chemical.
#'
#' @author Caroline Ring, John Wambaugh
#' @export fit_all
#' @importFrom PK nca.batch
#' @importFrom magrittr "%>%"

fit_all <- function(data.set,
                    model,
                    modelfun = NA,
                    ratio.data.to.dose = 1,

                    compound.col = "Compound",
                    dtxsid.col = "DTXSID",
                    cas.col = "CAS",

                    reference.col = "Reference",

                    species.col = "Species",
                    species.default = NULL,

                    species.weight.col = "Species.Weight",
                    species.weight.units.col = "Species.Weight.Units",
                    species.weight.units.default = NULL,

                    dose.col = "Dose",

                    time.col = "Time",
                    time.units.col = "Time.Units",
                    time.units.default = NULL,

                    media.col = "Media",
                    media.units.col = "Media.Units",
                    media.units.default = NULL,

                    value.col = "Value",

                    units.col = "Units",
                    units.default = NULL,

                    route.col = "Route",
                    route.default = NULL,

                    source.col = "Source",
                    source.default = NULL,

                    loq.col = "LOQ",
                    loq.default = NULL,

                    subject.col = "Subject",
                    subject.default = NULL,

                    info.col = "info",
                    info.default = NULL) {

  data.set <- rename_columns(data.set,
                             compound.col,
                             dtxsid.col,
                             cas.col,
                             reference.col,
                             species.col,
                             species.default,
                             species.weight.col,
                             species.weight.units.col,
                             species.weight.units.default,
                             dose.col,
                             time.col,
                             time.units.col,
                             time.units.default,
                             media.col,
                             media.units.col,
                             media.units.default,
                             value.col,
                             units.col,
                             units.default,
                             route.col,
                             route.default,
                             source.col,
                             source.default,
                             loq.col,
                             loq.default,
                             subject.col,
                             subject.default,
                             info.col,
                             info.default)

  ### Coerce all 'Value' values to be numeric
  data.set$Value <- as.numeric(data.set$Value)

  ### number of rows in data.set
  N.PREV <- dim(data.set)[1]

  ### display messages describing loaded data
  cat(paste(N.PREV, "concentration vs. time observations loaded.\n"))
  cat(paste(length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references.\n"))

  ### coerce 'Dose' values to numeric and say so
  if (is.character(class(data.set$Dose)))
  {
    cat("Column \"Dose\" converted to numeric.\n")
    data.set$Dose <- as.numeric(data.set$Dose)
  }

  ### coerce data.set to data.table object
  data.set <- as.data.table(data.set)

  # Right now code only recognizes "po" and "iv" as routes (my bad):
  ### coerce route names, 'oral' and 'intravenous', to 'po' and 'iv'
  data.set[Route == "oral", Route := "po"]
  data.set[Route == "intravenous", Route:="iv"]

  ### subset to data of routes 'po' and 'iv' only
  data.set <- data.set[Route %in% c("po", "iv")]
  cat(paste("Restricting to intravenous and oral routes eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))

  ### recalcuate N.PREV after removal of rows corresponding to other routes
  N.PREV <- dim(data.set)[1]

  # Harmonize the compound names:
  ### make all compound names completely lower-case
  data.set[, Compound := tolower(Compound)]

  ### do the same thing for species
  data.set[, Species := tolower(Species)]

  # This way the weight units cancel (must still pay attention to denominator
  # of data to determine units for Vd):

  ### coerce any 'Value' values of 0 to NA
  data.set$Value[data.set$Value == 0] <- NA
  cat("Converting 'Value' values of 0 to NA \n")

  ### normalize 'Value' by ratio.data.to.dose
  data.set[, Value := Value * ratio.data.to.dose]

  #Ignore data close to LOQ:
  data.set[Value < 2 * LOQ, Value := NA]

  #set an iv variable to TRUE/FALSE
  data.set[Route == 'iv', iv := TRUE]
  data.set[Route !='iv', iv := FALSE]

  #convert time from hours to days
  data.set <- data.set[!is.na(Time)]
  cat(paste("Requiring time to have a value != NA eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))
  N.PREV <- dim(data.set)[1]
  data.set[, Time.Days := Time / 24]
  data.set[, c('Max.Time.Days',
               'Time.Steps.PerHour') := list(max(Time.Days),
                                             1 / min(diff(c(0, sort(unique(Time)))))),
           by = .(DTXSID, Dose, Route)]

  # How many >LOQ observations do we have per chemical/species/reference?
  data.set[, N.Obs.Ref := dim(subset(.SD, !is.na(Value)))[1], by = .(Reference, DTXSID, Species)]
  # Not much we can do if fewer than 4 points (for instance, can't estimate Sigma'):
  data.set[, Usable := N.Obs.Ref > 3, by = .(DTXSID, Reference, Species, Route)]
  data.set <- data.set[Usable == TRUE]
  cat(paste("Restricting to references with more than three observations above LOQ eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))
  N.PREV <- dim(data.set)[1]

  data.set[,
           Usable := ifelse(Time >= .SD[Value == max(Value, na.rm = T), Time] |
                              !is.na(Value), T, F),
           by=.(Route,Reference,DTXSID,Species)]

  ### subset data to where Usable == TRUE
  data.set <- data.set[Usable == TRUE]
  cat(paste("Eliminating observations for doses that are below LOQ before the peak conc. is reached eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))

  ### recalcuate N.PREV after removal FALSE's
  N.PREV <- dim(data.set)[1]

  #############################
  #############################
  #############################
  # browser()
  #Non-comapartmental fits:
  if (model=="noncompartment") {

    PK.fit.table <- do_noncomp_fit(data.set)

  } else { ### but doesn't the else here imply that the model is not noncomp?
    ### should below comment say, "if model is not noncomp"?
    ### Okay, I'm pretty sure I'm correct, as params.by.cas.spec is called by 1-comp

    ### if model is noncomp, create data.frame of unique pairings of CAS, Species, and Compound
    params.by.cas.spec <- data.set[,unique(.SD[,.(Compound)]),by=.(DTXSID,CAS,Species)]

    #get the rat parameters for the 1-compartment model for each chemical
    if (model == '1compartment') {
      # browser()
      #      if (any(params.by.cas.spec$CAS%in%get_cheminfo(model='1compartment')))
      #        params.by.cas.spec[CAS %in% get_cheminfo(model='1compartment',species=Species),
      #        c("kelim", "Vdist", "Fgutabs", "kgutabs") :=
      #        httk::parameterize_1comp(chem.cas=CAS,
      #        default.to.human=TRUE,
      #        species=Species)[c("kelim","Vdist","Fgutabs","kgutabs")],
      #        by=c("CAS","Species")]

      # Wambaugh et al. (2018) medians:
      # apply(chem.invivo.PK.aggregate.data,2,function(x) median(as.numeric(x),na.rm=T))
      ### starting point for predictions
      params.by.cas.spec[, kelim := 0.25]
      params.by.cas.spec[, Vdist := 5.56]
      params.by.cas.spec[, Fgutabs := 1.0]
      params.by.cas.spec[, kgutabs := 2.19]

      #      if (modelfun=="analytic") params.by.cas.spec[, setdiff(names(params.by.cas.spec),
      #                                                        c("CAS",
      #                                                          "kelim",
      #                                                          "Vdist",
      #                                                          "Fgutabs",
      #                                                          "kgutabs")):=NULL]
    } else if (model == '2compartment') {
      # Use this when parameterize_2comp is implemented in httk
      #   params.by.cas <- data.set[, httk::parameterize_2comp(chem.cas=CAS,
      #                                                  default.to.human=TRUE,
      #                                                  species="Rat"),
      #                             by=CAS]
      if (any(params.by.cas.spec$DTXSID %in% httk::get_cheminfo(info="dtxsid")))
        params.by.cas.spec[DTXSID %in% httk::get_cheminfo(info="dtxsid"), ### where did get_cheminfo come from? httk?
                           c("kelim",
                             "Vdist",
                             "Fgutabs",
                             "kgutabs") := httk::parameterize_1comp(dtxsid = DTXSID,
                                                                    default.to.human = TRUE,
                                                                    species = ifelse(tolower(Species) %in% colnames(httk::physiology.data),
                                                                                     Species,
                                                                                     "Human"))[c("kelim",
                                                                                         "Vdist",
                                                                                         "Fgutabs",
                                                                                         "kgutabs")],
                           by = c("DTXSID", "Species")]
      # Wambaugh et al. (2018) medians:
      # apply(chem.invivo.PK.aggregate.data,2,function(x) median(as.numeric(x),na.rm=T))
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), kelim := 0.25]
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), Vdist := 5.56]
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), Fgutabs := 1.0]
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), kgutabs := 2.19]
      #      params.by.cas.spec <- tmp[, .(CAS, Vdist,kelim, Rblood2plasma,
      #                               MW, hematocrit, million.cells.per.gliver)]
      #Since we don't yet have 2comp predictions,
      #just choose some arbitrary numbers as starting points for the fit.
      params.by.cas.spec[, V1 := 1]
      params.by.cas.spec[, Ralphatokelim := 1.2]
      params.by.cas.spec[, Fbetaofalpha := 0.8]

      #      if (modelfun=="analytic") params.by.cas.spec[, setdiff(names(params.by.cas.spec),
      #                                                        c("CAS",
      #                                                          "kelim",
      #                                                          "Ralphatokelim",
      #                                                          "Fbetaofalpha",
      #                                                          "V1",
      #                                                          "Fgutabs",
      #                                                          "kgutabs")):=NULL]
    } else if (model == "flat") {

      params.by.cas.spec[, A := 1]
    }

    #############################
    #############################
    #############################

    ### merge data.set and params.by.cas.spec to add parameter columns
    data.set <- merge(data.set, params.by.cas.spec, by = c("Compound",
                                                           "DTXSID",
                                                           "CAS",
                                                           "Species"))
    ### create vector of parameter names using column names of params.by.cas.spec
    paramnames <- names(params.by.cas.spec)
    ### remove names that were only used for merging, i.e. CAS, Species, and Compound
    paramnames <- paramnames[!(paramnames %in% c("DTXSID",
                                                 "CAS",
                                                 "Species",
                                                 "Compound"))]

    #Replace spaces in references with "."
    data.set[, Reference := gsub(Reference,
                                 pattern = ' ',
                                 replacement = '.')]

    ### PK.fit.joint is a data.frame containing a row of parameter values per param.value.type per CAS
    PK.fit.joint <- data.set[,
                             analyze_pk_data(fitdata = .SD, ### what is .SD
                                             this.dtxsid = DTXSID,
                                             paramnames = paramnames,
                                             modelfun = modelfun,
                                             model = model),
                             by = c("DTXSID", "Species")]

    ### that correspond to multiple references
    multi.ref.cas <- data.set[, length(unique(Reference)) > 1, by = c("DTXSID", "Species")]
    multi.ref.cas <- subset(multi.ref.cas, V1 == T)$DTXSID
    if (length(multi.ref.cas) > 0) {
      data.set.multi.ref <- subset(data.set, DTXSID %in% multi.ref.cas)

      PK.fit.separate <- data.set.multi.ref[,
                                            analyze_pk_data(fitdata = .SD,
                                                            this.dtxsid = DTXSID,
                                                            paramnames = paramnames,
                                                            modelfun = modelfun,
                                                            model = model,
                                                            this.reference = Reference),
                                            by = c("DTXSID", "Species", "Reference")]

      #       PK.fit.separate.geomean <- PK.fit.separate %>% filter(param.value.type == "Fitted geometric mean")
      #       PK.fit.joint.geomean <- PK.fit.joint %>% filter(param.value.type == "Fitted geometric mean")
      #
      #       sigmas.joint <- PK.fit.joint.geomean[, grep("sigma", names(PK.fit.joint.geomean)), with = FALSE]
      #       sigmas.sep <- PK.fit.separate.geomean[, grep("sigma", names(PK.fit.separate.geomean)), with = FALSE]
      #
      #       sigmas.vec <- c(as.vector(sigmas.joint), as.vector(sigmas.sep))
      #
      #       PK.fit.joint <- PK.fit.joint[, -grep("sigma", names(PK.fit.joint)), with = FALSE]
      #       PK.fit.separate <- PK.fit.separate[, -grep("sigma", names(PK.fit.separate)), with = FALSE]

      PK.fit.bind <- rbind(PK.fit.joint,
                           PK.fit.separate, fill = TRUE)
    } else PK.fit.bind <- PK.fit.joint

    #record which model was fit to data and whether it was full or analytical
    PK.fit.bind[, model := model]
    PK.fit.bind[, model.type := modelfun]

    #############################
    #############################
    #############################

    #Do post fit calculations:
    ### basically calculate parameters based on previously calculated paramaters
    PK.fit.table <- PK.fit.bind
    if (model == "1compartment") {
      #Get stats for fitted total clearance :    L/kg body weight/h
      PK.fit.table[, CLtot := Vdist * kelim]

      #      browser()
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
                   halflife := log(2) / kelim]


      PK.fit.table[,
                   tpeak.oral := log(kgutabs / kelim) / (kgutabs - kelim)]

      ### why does Fgutabs have to be 1 right here
      ### this just copies the Ccompartment column
      ### not sure how values correspond to specific param.value.types
      ### again, it's just a time, conc, and auc matrix
      PK.fit.table[,
                   Cpeak.oral.1mgkg := analytic_1comp_fun(
                     params=list(Fgutabs = 1,
                                 kgutabs = kgutabs,
                                 kelim = kelim,
                                 Vdist = Vdist),
                     dose = 1,
                     tpeak.oral,
                     iv.dose = F)[, "Ccompartment"]]

      ### subset only certain param.value.types, not sure why
      PK.fit.table <- PK.fit.table[param.value.type %in% c("Predicted",
                                                           "Fitted geometric mean",
                                                           "Fitted geometric std dev")]

    } else if(model == "2compartment") {
      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   c("beta",
                     "alpha") := list(Fbetaofalpha*Ralphatokelim*kelim,
                                      Ralphatokelim*(kelim+10^-6))]

      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   c("k21",
                     "k12") := list(alpha * beta / kelim,
                                    alpha + beta - kelim - alpha * beta / kelim)]

      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   c("halflife",
                     "Vss",
                     "CLtot",
                     "Varea.or.Vbeta") := list(log(2) / beta,
                                               V1 * (k21 + k12) / k21,
                                               V1 * (k21 + k12) / k21 * beta,
                                               V1 * kelim / beta)]
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
      PK.fit.table[param.value.type != "Predicted",
                   Css := (Fgutabs * 1) / (CLtot * 24)]

      PK.fit.table <- PK.fit.table[param.value.type %in% c("Predicted",
                                                           "Fitted geometric mean",
                                                           "Fitted geometric std dev")]
    } else if (model == 'flat') {
      PK.fit.table <- PK.fit.table[param.value.type %in% c("Predicted",
                                                           "Fitted geometric mean",
                                                           "Fitted geometric std dev")]
    }
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



