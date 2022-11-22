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
#' @param LOQ_factor Numeric. Observations with concentrations less than
#'   `LOQ_factor * LOQ` will be removed. Default 2.
#' @param get_starts_args Any additional arguments to [get_starts()] (other than
#'   `model` and `fitdata`, which are always passed). Default NULL to accept the
#'   default arguments for [get_starts()].
#' @param get_lower_args Any additional arguments to [get_lower_bounds()] (other
#'   than `model` and `fitdata`, which are always passed). Default NULL to
#'   accept the default arguments for [get_lower_bounds()].
#' @param get_upper_args Any additional arguments to [get_upper_bounds()] (other
#'   than `model` and `fitdata`, which are always passed). Default NULL to
#'   accept the default arguments for [get_upper_bounds()].
#' @param suppress.messages Logical: Whether to suppress verbose messages.
#'   Default FALSE, to be verbose.
#'
#' @return A `data.table` of fitted parameter values for each chemical.
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

                    LOQ_factor = 2,

                    get_starts_args = NULL,
                    get_lower_args = NULL,
                    get_upper_args = NULL,
                    optimx_args = list(
                      "method" = "L-BFGS-B",
                      "control" = list("factr" = 1e7,
                                       "fnscale" = -1)
                    ),

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
                             info.default,

                             LOQ_factor = LOQ_factor,
                             suppress.messages = suppress.messages)

  ### recalcuate N.PREV after preprocessing
  N.PREV <- dim(data.set)[1]


  #Non-comapartmental fits:
  if (model=="noncompartment") {

    PK.fit.table <- do_noncomp_fit(data.set)

  } else {

    #Analyze by chemical & species first.
if(!suppress.messages){
  message("Analyzing data by chemical and species, with each reference having its own error SD...")
}

    # tmp <- split(data.set,
    #             by = c("DTXSID", "Species"))
    # #
    # PK_fit_joint_list <- lapply(tmp,
    #                            function(fitdata){
    #                              analyze_subset(fitdata = fitdata,
    #                                             #fitdata = .SD passes a data.table
    #                                             #consisting of the columns
    #                                             #referenced in .SDcols, split by the
    #                                             #values in "by"
    #                                             modelfun = modelfun,
    #                                             model = model,
    #                                             pool_sigma = FALSE,
    #                                             LOQ_factor = LOQ_factor,
    #                                             get_starts_args = get_starts_args,
    #                                             get_lower_args = get_lower_args,
    #                                             get_upper_args = get_upper_args,
    #                                             optimx_args = optimx_args,
    #                                             suppress.messages=suppress.messages)
    #                            } )

    PK.fit.joint <- data.set[,

                             analyze_subset(fitdata = .SD,
                                            #fitdata = .SD passes a data.table
                                            #consisting of the columns
                                            #referenced in .SDcols, split by the
                                            #values in "by"
                                             modelfun = modelfun,
                                             model = model,
                                            pool_sigma = FALSE,
                                            LOQ_factor = LOQ_factor,
                                            get_starts_args = get_starts_args,
                                            get_lower_args = get_lower_args,
                                            get_upper_args = get_upper_args,
                                            optimx_args = optimx_args,
                                             suppress.messages=suppress.messages),
                             .SDcols = names(data.set),
                             #.SDcols argument: By default, .SD includes all
                             #columns *except* those named in "by=". We want to
                             #*include* those named in "by=", so we explicitly
                             #state that .SDcols must include all columns in
                             #data.set. Otherwise, "fitdata" would be missing
                             #the "DTXSID" and "Species" columns.
                             by = c("DTXSID", #split data into unique combinations of these columns
                                    "Species")
                             ]

    ### If chemical/species combination has multiple references:
    ### then analyze each reference individually,
    ###and also do a pooled analysis without individual reference sigmas
    data.set[, MultipleReferences := length(unique(Reference)) > 1,
             by = c("DTXSID", "Species")]
    multi.ref.cas <- unique(subset(data.set,
                                   MultipleReferences == TRUE)$DTXSID)
    if (length(multi.ref.cas) > 0) {
      if(!suppress.messages){
        message("Analyzing data by chemical, species, and reference...")
      }
      data.set.multi.ref <- subset(data.set, MultipleReferences == TRUE)
      data.set.multi.ref[, Reference_orig := Reference]
      PK.fit.separate <- data.set.multi.ref[,
                                 analyze_subset(fitdata = .SD,
                                                 modelfun = modelfun,
                                                 model = model,
                                                pool_sigma = FALSE,
                                                LOQ_factor = LOQ_factor,
                                                get_starts_args = get_starts_args,
                                                get_lower_args = get_lower_args,
                                                get_upper_args = get_upper_args,
                                                optimx_args = optimx_args,
                                                suppress.messages = suppress.messages),
                                 .SDcols = names(data.set.multi.ref),
                                                 by = c(
                                                   "DTXSID",
                                                   "Species",
                                                   "Reference_orig")]


      if(!suppress.messages){
        message("Analyzing data by chemical and species, pooling all references...")
      }
      PK.fit.pooled <- data.set.multi.ref[,
                                          analyze_subset(fitdata = .SD,
                                                         modelfun = modelfun,
                                                         model = model,
                                                         pool_sigma = TRUE,
                                                         LOQ_factor = LOQ_factor,
                                                         get_starts_args = get_starts_args,
                                                         get_lower_args = get_lower_args,
                                                         get_upper_args = get_upper_args,
                                                         optimx_args = optimx_args,
                                                         suppress.messages = suppress.messages),
                                          .SDcols = names(data.set.multi.ref),
                                          by = c(
                                            "DTXSID",
                                            "Species")]

      PK.fit.bind <- rbindlist(list("Joint" = PK.fit.joint,
                           "Separate" = PK.fit.separate,
                            "Pooled" = PK.fit.pooled
      ),
                           use.names = TRUE,
                           fill = TRUE,
                           idcol = "Analysis_Type")
    } else PK.fit.bind <- PK.fit.joint

    #record which model was fit to data and whether it was full or analytical
    PK.fit.bind[, model := model]
    PK.fit.bind[, model.type := modelfun]

    #############################
    #############################
    #############################
}

  return(PK.fit.bind)
}
