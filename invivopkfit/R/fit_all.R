#' Main fitting function
#'
#' Fits parameters of a specified model to concentration-time data given in
#' data.set
#'
#' @param data.set A \code{data.frame} of concentration-time data. Preferably
#'   \code{pkdataset_nheerlcleaned} or \code{pkdataset_nheerlorig}.
#' @param model Which general model should be fit for each chemical. Presently,
#'   only "1compartment" and "2compartment" are implemented.
#' @param modelfun Either "analytic" or "full" -- whether to fit using the
#'   analytic solution to the model, or the full ODE model. Presently,
#'   "analytic" is recommended (because the analytic solution is exact and much
#'   faster).
#' @param ratio.data.to.dose Ratio between the mass units used to report the
#'   concentration data and the mass units used to report the dose. Default 1.
#'   For example, concentration reported in ug/L and dose reported in mg/kg/day
#'   would require \code{ratio.data.to.dose = 0.001}, because 1 ug/1 mg = 1e-6 g
#'   / 1e-3 g = 0.001.
#' @param compound.col Column name in \code{data.set} that identifies chemical
#'   compound. Default "Compound".
#' @param cas.col Column name in \code{data.set} to identify CASRN.
#'   Default"CAS".
#' @param reference.col Column name in \code{data.set} to identify reference.
#'   Default "Reference".
#' @param species.col Column name in \code{data.set} to identify species.
#'   Default "Species".
#' @param species.default If no species column exists in \code{data.set}, one
#'   will be created and filled with this value.  Default NULL.
#' @param species.weight.col Column name in \code{data.set} to identify species
#'   weight. Default "Species.Weight".
#' @param species.weight.units.col Column name in \code{data.set} to identify
#'   species weigh units. Default"Species.Weight.Units".
#' @param species.weight.units.default If no species weight units column exists
#'   in \code{data.set}, one will be created and filled with this value. Default
#'   NULL.
#' @param dose.col Column name in \code{data.set} to identify dose. Default
#'   "Dose".
#' @param time.col Column name in \code{data.set} to identify time. Default
#'   "Time".
#' @param time.units.col Column name in \code{data.set} to identify time units.
#'   Default "Time.Units."
#' @param time.units.default If no time units column exists in \code{data.set},
#'   one will be created and filled with this value.  Default NULL.
#' @param media.col Column name in \code{data.set} to identify media. Default
#'   "Media".
#' @param media.units.col Column name in \code{data.set} to identify media
#'   units. Default "Media.Units".
#' @param media.units.default If no media units column exists in
#'   \code{data.set}, one will be created and filled with this value.  Default
#'   NULL.
#' @param value.col Column name in \code{data.set} to identify value. Default
#'   "Value".
#' @param units.col Column name in \code{data.set} to identify units. Default
#'   "Units".
#' @param units.default If no units column exists in \code{data.set}, one will
#'   be created and filled with this value.  Default NULL.
#' @param route.col Column name in \code{data.set} to identify route of
#'   administration. Default "Route".
#' @param route.default If no route column exists in \code{data.set}, one will
#'   be created and filled with this value.  Default NULL.
#' @param source.col Column name in \code{data.set} to identify source. Default
#'   "Source."
#' @param source.default If no source column exists in \code{data.set}, one will
#'   be created and filled with this value.  Default NULL.
#' @param loq.col Column name in \code{data.set} to identify LOQ. Default "LOQ".
#' @param loq.default If no LOQ column exists in \code{data.set}, one will be
#'   created and filled with this value.  Default NULL.
#' @param subject.col Column name in \code{data.set} to identify subject.
#'   Default "Subject."
#' @param subject.default If no subject column exists in \code{data.set}, one
#'   will be created and filled with this value.  Default NULL.
#' @param info.col Column name in \code{data.set} to serve as info column.
#'   Default "Info".
#' @param info.default If no info column exists in \code{data.set}, one will be
#'   created and filled with this value.  Default NULL.
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
                    info.default = NULL,
                    suppress.messages = FALSE,
                    sig.figs=4) {

  #############################
  #############################
  #############################
  ## PREPROCESS DATA #########
  #############################
  #############################
  #############################

 data.set <- preprocess_data(data.set,
                             ratio.data.to.dose,

                             compound.col ,
                             dtxsid.col,
                             cas.col,

                             reference.col,

                             species.col ,
                             species.default,

                             species.weight.col ,
                             species.weight.units.col ,
                             species.weight.units.default ,

                             dose.col,

                             time.col,
                             time.units.col ,
                             time.units.default ,

                             media.col,
                             media.units.col ,
                             media.units.default,

                             value.col ,

                             units.col ,
                             units.default,

                             route.col,
                             route.default ,

                             source.col ,
                             source.default ,

                             loq.col ,
                             loq.default ,

                             subject.col ,
                             subject.default,

                             info.col,
                             info.default)

  ### recalcuate N.PREV after preprocessing
  N.PREV <- dim(data.set)[1]

  #############################
  #############################
  #############################
  ## SET STARTING GUESSES ####
  #############################
  #############################
  #############################

  #Non-comapartmental fits:
  if (model=="noncompartment") {

    PK.fit.table <- do_noncomp_fit(data.set)

  } else {
    ### if model is not noncomp, create data.frame of unique pairings of CAS, Species, and Compound
    params.by.cas.spec <- data.set[,
                                   unique(.SD[,.(Compound)]),
                                   by=.(DTXSID,CAS,Species)]
    if (model == '1compartment') {
      ### starting point for optimizer
      params.by.cas.spec[, kelim := 0.25]
      params.by.cas.spec[, Vdist := 5.56]
      params.by.cas.spec[, Fgutabs := 1.0]
      params.by.cas.spec[, kgutabs := 2.19]

    } else if (model == '2compartment') {
      if (any(params.by.cas.spec$DTXSID %in% httk::get_cheminfo(info="dtxsid"))){
        #if we have PK params in httk, use those as starting points.
        params.by.cas.spec[DTXSID %in% httk::get_cheminfo(info="dtxsid"),
                           c("kelim",
                             "Vdist",
                             "Fgutabs",
                             "kgutabs") := httk::parameterize_1comp(dtxsid = DTXSID,
                                                                    default.to.human = TRUE,
                                                                    suppress.messages = TRUE,
                                                                    species = ifelse(tolower(Species) %in% colnames(httk::physiology.data),
                                                                                     Species,
                                                                                     "Human"))[c("kelim",
                                                                                         "Vdist",
                                                                                         "Fgutabs",
                                                                                         "kgutabs")],
                           by = c("DTXSID", "Species")]
      }
     #if we don't have PK params in httk, set default starting points.
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), kelim := 0.25]
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), Vdist := 5.56]
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), Fgutabs := 1.0]
      params.by.cas.spec[!params.by.cas.spec$DTXSID %in%
                         httk::get_cheminfo(info="dtxsid"), kgutabs := 2.19]

      #Since we don't yet have 2comp predictions,
      #just choose some arbitrary numbers as starting points for the fit.
      params.by.cas.spec[, V1 := 1]
      params.by.cas.spec[, Ralphatokelim := 1.2]
      params.by.cas.spec[, Fbetaofalpha := 0.8]

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
    ### Throughout the code we make use of the data.table feature ".SD".
    ### .SD stands for "subset of the data". In a call to a data.table object A:
    ### A[i,j, by = k]
    ### i indicates the rows that are impacted (fuzzy on this) (that's correct -- CLR)
    ### j indicates a column that is changed
    ### k indicates
    ###
    ### For a complete explanation go to:
    ### https://cran.r-project.org/web/packages/data.table/vignettes/datatable-sd-usage.html
                             analyze_pk_data(fitdata = .SD,
                                             this.dtxsid = DTXSID,
                                             paramnames = paramnames,
                                             modelfun = modelfun,
                                             model = model,
                                             suppress.messages=suppress.messages),
                             by = c("DTXSID", "Species")]

    #browser()
    ### Rerun subsetting per reference just for chemical/species comvbinations
    ### that have multiple references:
    data.set[, MultipleReferences := length(unique(Reference)) > 1,
             by = c("DTXSID", "Species")]
    multi.ref.cas <- unique(subset(data.set,
                                   MultipleReferences == TRUE)$DTXSID)
    if (length(multi.ref.cas) > 0) {
      data.set.multi.ref <- subset(data.set, MultipleReferences == TRUE)

      PK.fit.separate <- data.set.multi.ref[,
                                 analyze_pk_data(fitdata = .SD,
                                                 this.dtxsid = DTXSID,
                                                 paramnames = paramnames,
                                                 modelfun = modelfun,
                                                 model = model,
                                                 this.reference = Reference,
                                                 suppress.messages = suppress.messages),
                                                 by = c(
                                                   "DTXSID",
                                                   "Species",
                                                   "Reference")]

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
      PK.fit.table[, CLtot := signif(Vdist * kelim, sig.figs)]

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



