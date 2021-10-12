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
#' @param compound.col placeholder
#' @param cas.col placeholder
#' @param reference.col placeholder
#' @param species.col placeholder
#' @param species.default placeholder
#' @param species.weight.col placeholder
#' @param species.weight.units.col placeholder
#' @param species.weight.units.default placeholder
#' @param dose.col placeholder
#' @param time.col placeholder
#' @param time.units.col placeholder
#' @param time.units.default placeholder
#' @param media.col placeholder
#' @param media.units.col placeholder
#' @param media.units.default placeholder
#' @param value.col placeholder
#' @param units.col placeholder
#' @param units.default placeholder
#' @param route.col placeholder
#' @param route.default placeholder
#' @param source.col placeholder
#' @param source.default placeholder
#' @param loq.col placeholder
#' @param loq.default placeholder
#' @param subject.col placeholder
#' @param subject.default placeholder
#' @param info.col placeholder
#' @param info.default placeholder
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

                    series.col = "Series",
                    series.default = NULL)
{

  if(!is.null(species.default)) {
    data.set[, species.col] <- species.default
  }

  if(!is.null(species.weight.units.default)) {
    data.set[, species.weight.units.col] <- species.weight.units.default
  }

  if(!is.null(time.units.default)) {
    data.set[, time.units.col] <- time.units.default
  }

  if(!is.null(media.units.default)) {
    data.set[, media.units.col] <- media.units.default
  }

  if(!is.null(units.default)) {
    data.set[, units.col] <- units.default
  }

  if(!is.null(route.default)) {
    data.set[, route.col] <- route.default
  }

  if(!is.null(loq.default)) {
    data.set[, loq.col] <- loq.default
  }

  if(!is.null(subject.default)) {
    data.set[, subject.col] <- subject.default
  }

  if(!is.null(info.default)) {
    data.set[, info.col] <- info.default
  }

  if(!is.null(series.default)) {
    data.set[, series.col] <- series.default
  }

  cols <- c(
    compound.col,
    cas.col,
    reference.col,
    species.col,
    species.weight.col,
    species.weight.units.col,
    dose.col,
    time.col,
    time.units.col,
    media.col,
    media.units.col,
    value.col,
    units.col,
    route.col,
    source.col,
    loq.col,
    subject.col,
    info.col,
    series.col
  )

  ### stop if any columns missing from data.set
  if (!(all(cols %in% colnames(data.set))))
  {
    stop(paste("Missing columns named:",
               paste(cols[!(cols%in%colnames(data.set))],collapse=", ")))
  }

  ### Set column order of data.table to cols vector
  data.set <- data.table::setcolorder(data.set, cols)

  # Standardize the column names:
  compound.col <- "Compound"
  cas.col <- "CAS"
  reference.col <- "Reference"
  species.col <- "Species"
  species.weight.col <- "Species.Weight"
  species.weight.units.col <- "Species.Weight.Units"
  dose.col <- "Dose"
  time.col <- "Time"
  time.units.col <- "Time.Units"
  media.col <- "Media"
  media.units.col <- "Media.Units"
  value.col <- "Value"
  units.col <- "Units"
  route.col <- "Route"
  source.col <- "Source"
  loq.col <- "LOQ"
  subject.col <- "Subject"
  info.col <- "info"
  series.col <- "Series.ID"

  ### rename colnames to standardized colnames
  colnames(data.set) <- c(
    compound.col,
    cas.col,
    reference.col,
    species.col,
    species.weight.col,
    species.weight.units.col,
    dose.col,
    time.col,
    time.units.col,
    media.col,
    media.units.col,
    value.col,
    units.col,
    route.col,
    source.col,
    loq.col,
    subject.col,
    info.col,
    series.col
  )

  if(any(data.set$Series.ID == "hack")) {
    setnames(data.set, "Series.ID", "delete")
    data.set[, "delete" := NULL]
  }

  ### rename extra/NA columns in template as "delete"
  names_corrected <- tidyr::replace_na(colnames(data.set), "delete")
  colnames(data.set) <- names_corrected

  ##############
  # make addition of series.col argument compatible with old data.set, chris. you have work to do here.
  ### delete columns named "delete"
  if("delete" %in% colnames(data.set)) {data.set[, "delete" == names(data.set)] <- NULL}

  ### Coerce all 'Value' values to be numeric
  data.set$Value <- as.numeric(data.set$Value)

  ### number of rows in data.set
  N.PREV <- dim(data.set)[1]

  ### display messages describing loaded data
  cat(paste(N.PREV,"concentration vs. time observations loaded.\n"))
  cat(paste(length(unique(data.set$CAS)),"unique chemicals,",
            length(unique(data.set$Species)),"unique species, and",
            length(unique(data.set$Reference)),"unique references.\n"))

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
  data.set[Route=="oral",Route:="po"]
  data.set[Route=="intravenous",Route:="iv"]

  ### subset to data of routes 'po' and 'iv' only
  data.set <- data.set[Route %in% c("po","iv")]
  cat(paste("Restricting to intravenous and oral routes eliminates",
            N.PREV - dim(data.set)[1],"observations.\n"))
  cat(paste(dim(data.set)[1],"observations of",
            length(unique(data.set$CAS)),"unique chemicals,",
            length(unique(data.set$Species)),"unique species, and",
            length(unique(data.set$Reference)),"unique references remain.\n"))

  ### recalcuate N.PREV after removal of rows corresponding to other routes
  N.PREV <- dim(data.set)[1]

  # Harmonize the compound names:
  ### make all compound names completely lower-case
  data.set[,Compound:=tolower(Compound)]

  ### do the same thing for species
  data.set[, Species := tolower(Species)]

  # This way the weight units cancel (must still pay attention to denominator
  # of data to determine units for Vd):

  ### coerce any 'Value' values of 0 to NA
  data.set$Value[data.set$Value == 0] <- NA
  cat("Converting 'Value' values of 0 to NA \n")

  ### normalize 'Value' by ratio.data.to.dose
  data.set[,Value:=Value*ratio.data.to.dose]

  #Ignore data close to LOQ:
  data.set[Value<2*LOQ,Value:=NA]

  #set an iv variable to TRUE/FALSE
  data.set[Route=='iv', iv:=TRUE]
  data.set[Route!='iv', iv:=FALSE]

  #convert time from hours to days
  data.set <- data.set[!is.na(Time)]
  cat(paste("Requiring time to have a value != NA eliminates",
            N.PREV - dim(data.set)[1],"observations.\n"))
  cat(paste(dim(data.set)[1],"observations of",
            length(unique(data.set$CAS)),"unique chemicals,",
            length(unique(data.set$Species)),"unique species, and",
            length(unique(data.set$Reference)),"unique references remain.\n"))
  N.PREV <- dim(data.set)[1]
  data.set[, Time.Days:=Time/24]
  data.set[, c('Max.Time.Days',
               'Time.Steps.PerHour'):=list(max(Time.Days),
                                           1/min(diff(c(0,sort(unique(Time)))))),
           by=.(CAS, Dose, Route)]

  # How many >LOQ observations do we have per chemical/species/reference?
  data.set[,N.Obs.Ref:=dim(subset(.SD,!is.na(Value)))[1],by=.(Reference,CAS,Species)]
  # Not much we can do if fewer than 4 points (for instance, can't estimate Sigma'):
  data.set[,Usable:=N.Obs.Ref>3,by=.(CAS,Reference,Species,Route)]
  data.set <- data.set[Usable==TRUE]
  cat(paste("Restricting to references with more than three observations above LOQ eliminates",
            N.PREV - dim(data.set)[1],"observations.\n"))
  cat(paste(dim(data.set)[1],"observations of",
            length(unique(data.set$CAS)),"unique chemicals,",
            length(unique(data.set$Species)),"unique species, and",
            length(unique(data.set$Reference)),"unique references remain.\n"))
  N.PREV <- dim(data.set)[1]

  data.set[,
           Usable:=ifelse(Time >= .SD[Value==max(Value,na.rm=T),Time] |
                            !is.na(Value),T,F),
           by=.(Route,Reference,CAS,Species)]

  ### subset data to where Usable == TRUE
  data.set <- data.set[Usable==TRUE]
  cat(paste("Eliminating observations for doses that are below LOQ before the peak conc. is reached eliminates",
            N.PREV - dim(data.set)[1],"observations.\n"))
  cat(paste(dim(data.set)[1],"observations of",
            length(unique(data.set$CAS)),"unique chemicals,",
            length(unique(data.set$Species)),"unique species, and",
            length(unique(data.set$Reference)),"unique references remain.\n"))

  ### recalcuate N.PREV after removal FALSE's
  N.PREV <- dim(data.set)[1]

  #############################
  #############################
  #############################

  ### function to take average of multiple timepoints
  avg_value_fun <- function(data) {

    ### convert to data.frame, but fix this later because converting back and forth is silly
    data <- as.data.frame(data)

    ### if there are replicate timepoints, take the average Value for each timepoint
    if(length(unique(data$Time)) < length(data$Time)) {
      data <- data %>%
        group_by(Time) %>%
        mutate(mean = mean(Value, na.rm = TRUE))

      ### if all values are NA, output will return 'NaN's
      ### coerces 'NaN's back to NA
      data$mean[is.nan(data$mean)] <- NA

      ### remove replicate timepoints
      data <- data %>%
        distinct(Time, .keep_all = TRUE)

      ### delete Value column and reassign new mean column as Value column
      data$Value <- NULL
      names(data)[names(data) == "mean"] <- "Value"
    }

    data <- as.data.table(data)

    data
  }

  data.list <- list()
  for(this.compound in unique(data.set[, Compound])) {
    this.compound.data <- data.set[Compound == this.compound]

    ### list of columns by which to split data.set
    col_list <- list(this.compound.data$Reference,
                       this.compound.data$Dose,
                       this.compound.data$Route,
                       this.compound.data$Media)

    ### split data.set into list of data.tables
    split_list <- split(this.compound.data, col_list)

    ### remove data.tables with 0 rows
    # this.medium.list <- keep(this.medium.list, ~nrow(.) > 0)

    split_list_avg <- lapply(split_list, avg_value_fun)
    test <- do.call(rbind, split_list_avg)

    data.list[[this.compound]] <- test
  }

  ### combine data.tables back into one large data.set
  data.set <- do.call(rbind, data.list)

  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################

  #Non-comapartmental fits:
  if (model=="noncompartment")
  {
    ### average conc over multiple replicates/timepoints by study, dose, compound
    ### (in cvtdb study incorporates all said)

    # Normalize all concentrations by dose to fit jointly:
    data.set[,Value.Norm:=Value/Dose]

    ### not sure why PK.fit.table needs to be defined as NULL right here
    PK.fit.table<- NULL

    ### this is where the data is waded through and this.route.subset is defined
    for (this.cas in sort(unique(data.set$CAS)))
    {
      this.subset <- subset(data.set,CAS==this.cas&!is.na(Value.Norm))
      for (this.species in sort(unique(this.subset$Species)))
      {
        # Need oral and iv to estimate bioavailability:
        if ("po" %in% unique(this.subset$Route) & "iv" %in% unique(this.subset$Route))
        {
          this.iv.time.list <- list()
          this.iv.conc.list <- list()
          this.iv.id.list <- list()
          this.oral.time.list <- list()
          this.oral.conc.list <- list()
          this.oral.id.list <- list()
          this.row <- data.frame(Compound=this.subset$Compound[1],
                                 CAS=this.cas,
                                 Reference=NA,
                                 Source=NA,
                                 Species=this.species,
                                 Mean.iv.Dose=mean(subset(this.subset,Route=="iv")$Dose), ### averaged doses for iv and po
                                 Mean.po.Dose=mean(subset(this.subset,Route=="po")$Dose),
                                 model.type="noncompartmental",
                                 param.value.type="Estimated",
                                 AUC.po=NA,
                                 Vd.po=NA,
                                 CL.po=NA,
                                 AUC.iv=NA,
                                 Vd.iv=NA,
                                 CL.iv=NA,
                                 Fbio=NA,
                                 stringsAsFactors=F)

          this.row <- rbind(this.row,this.row)

          rownames(this.row) <- NULL
          this.row[2,"param.value.type"] <- "Estimated std dev"
          for (this.route in unique(this.subset$Route))
          {
            this.route.subset <- subset(this.subset,Route==this.route)

            new.route.subset <- NULL

            if (any(is.na(this.route.subset$Subject)))
            {
              for (this.reference in unique(this.route.subset$Reference))
              {
                for (this.dose in unique(this.route.subset$Dose))
                {
                  this.route.subset[is.na(this.route.subset$Subject)&this.route.subset$Dose==this.dose,"Subject"] <- paste(this.reference,this.route,this.dose,sep="-")
                  this.dose.subset <- subset(this.route.subset, Dose == this.dose)

                  ### make Value.Norm the average of all Value.Norm values
                  ### nest replicate data (Value) and have each row be a unique Time and Value.Norm
                  browser()
                  this.dose.subset <- this.dose.subset %>%
                    dplyr::group_by(Time) %>%
                    dplyr::mutate(Value.Norm = mean(Value.Norm)) %>%
                    dplyr::group_by_at(setdiff(names(this.route.subset), "Value")) %>%
                    tidyr::nest()

                  new.route.subset <- rbind(new.route.subset, this.dose.subset)
                }
              }
            }
            ### test this
            this.route.subset <- new.route.subset

            ### conc and time are required field for the PK package, which the non-compartmental model relies on
            this.route.subset$id <- paste(this.row[1,"Reference"],this.route.subset$Subject,sep="-")
            this.route.subset$conc <- this.route.subset$Value.Norm
            this.route.subset$time <- this.route.subset$Time
            for (this.id in unique(this.route.subset$id))
            {
              ### what does this do?
              if (dim(subset(this.route.subset,id==this.id))[1]<3)
              {
                this.route.subset <- subset(this.route.subset,id!=this.id)
              }
            }

            ### order this.route.subset by time
            this.route.subset <- this.route.subset[order(this.route.subset$time),]

            # if there are multiple series in a study, average concs at each time point BEFORE fit_all
            # i don't think the previous comment is correct; as the original data.set appears to have replicates

            # remove rows where 'time' value is not in the group of max number of unique time values
            # need to go back and see if original data.set has data that does not abide by this rule of having equal lengths of unique values
            # browser()
            # if (length(unique(table(this.route.subset$time)))>1) {
            #   this.route.subset.gdf <- this.route.subset %>% group_by(time) %>% tally()
            #   this.route.subset <- full_join(this.route.subset, this.route.subset.gdf)
            #   this.route.subset <- this.route.subset[n == max(n)]
            #   this.route.subset <- this.route.subset[,n:=NULL]
            # }

            if (dim(this.route.subset)[1]>2)
            {
              if (this.route=="po")
              {
                dose.arg <- 0
              }
              else if (this.route=="iv")
              {
                dose.arg <- 1
              }

              ### this tells which type of nca function to try, .complete or .batch
              if (length(this.route.subset$time)==length(unique(this.route.subset$time)))
              {
                out <- try(PK::nca.complete(data=this.route.subset,dose=dose.arg,method="z"))
              }
              else
              {
                out <- try(PK::nca.batch(data=this.route.subset,dose=dose.arg,method="z"))
              }

              ### not sure what's happening here, but if the code hits browser, something went wrong
              if (class(out)!="try-error")
              {
                out <- out$CIs
                if (!is.na(out[1,"est"]>0)) if (out[1,"est"]>0)
                {
                  this.row[1,paste("AUC",this.route,sep=".")] <-out[1,"est"]
                  this.row[2,paste("AUC",this.route,sep=".")] <-out[1,"stderr"]
                }
                if (!is.na(out[6,"est"]>0)) if (out[6,"est"]>0)
                {
                  this.row[1,paste("CL",this.route,sep=".")] <-out[6,"est"]
                  this.row[2,paste("CL",this.route,sep=".")] <-out[6,"stderr"]
                }
                if (!is.na(out[7,"est"]>0)) if (out[7,"est"]>0)
                {
                  this.row[1,paste("Vd",this.route,sep=".")] <-out[7,"est"]
                  this.row[2,paste("Vd",this.route,sep=".")] <-out[7,"stderr"]
                }
              } else browser()
            }
          }
          if (!is.na(this.row[1,"AUC.po"]) & !is.na(this.row[1,"AUC.iv"]))
          {
            this.row[1,"Fbio"] <- this.row[1,"AUC.po"]/this.row[1,"AUC.iv"]
            this.row[2,"Fbio"] <- this.row[1,"Fbio"]*((this.row[2,"AUC.po"]/this.row[1,"AUC.po"])^2+(this.row[2,"AUC.iv"]/this.row[1,"AUC.iv"])^2)^(1/2)
          }

          ### if PK.fit.table is null until here, why define it earlier?
          PK.fit.table <- rbind(PK.fit.table,this.row)
        }
      }
    }

    ### change PK.fit.table to data.table
    PK.fit.table <- as.data.table(PK.fit.table)

    ### check this Chris. Are changes in dose level accounted for?
    ### Current sets have multiple dose levels/set of concentration~time.

    # Multiply by dose to return from normalized uinits:
    PK.fit.table[,Vd.iv:=Vd.iv*Mean.iv.Dose]
    PK.fit.table[,Vd.po:=Vd.po*Mean.po.Dose]

  } else { ### but doesn't the else here imply that the model is not noncomp?
    ### should below comment say, "if model is not noncomp"?
    ### Okay, I'm pretty sure I'm correct, as params.by.cas.spec is called by 1-comp

    ### if model is noncomp, create data.frame of unique pairings of CAS, Species, and Compound
    params.by.cas.spec <- data.set[,unique(.SD[,.(Compound)]),by=.(CAS,Species)]

    ###########################################################################################################
    ###########################################################################################################
    ###########################################################################################################

    #get the rat parameters for the 1-compartment model for each chemical
    if (model=='1compartment')
    {
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
      params.by.cas.spec[, kelim:=0.25]
      params.by.cas.spec[, Vdist:=5.56]
      params.by.cas.spec[, Fgutabs:=1.0]
      params.by.cas.spec[, kgutabs:=2.19]

      #      if (modelfun=="analytic") params.by.cas.spec[, setdiff(names(params.by.cas.spec),
      #                                                        c("CAS",
      #                                                          "kelim",
      #                                                          "Vdist",
      #                                                          "Fgutabs",
      #                                                          "kgutabs")):=NULL]
    } else if (model=='2compartment') {
      # Use this when parameterize_2comp is implemented in httk
      #   params.by.cas <- data.set[, httk::parameterize_2comp(chem.cas=CAS,
      #                                                  default.to.human=TRUE,
      #                                                  species="Rat"),
      #                             by=CAS]
      if (any(params.by.cas.spec$CAS %in% get_cheminfo()))
        params.by.cas.spec[CAS %in% get_cheminfo(), ### where did get_cheminfo come from? httk?
                           c("kelim", "Vdist", "Fgutabs", "kgutabs") :=
                             httk::parameterize_1comp(chem.cas=CAS,
                                                      default.to.human=TRUE,
                                                      species=Species)[c("kelim","Vdist","Fgutabs","kgutabs")],
                           by=c("CAS","Species")]
      # Wambaugh et al. (2018) medians:
      # apply(chem.invivo.PK.aggregate.data,2,function(x) median(as.numeric(x),na.rm=T))
      params.by.cas.spec[!params.by.cas.spec$CAS%in%get_cheminfo(), kelim:=0.25]
      params.by.cas.spec[!params.by.cas.spec$CAS%in%get_cheminfo(), Vdist:=5.56]
      params.by.cas.spec[!params.by.cas.spec$CAS%in%get_cheminfo(), Fgutabs:=1.0]
      params.by.cas.spec[!params.by.cas.spec$CAS%in%get_cheminfo(), kgutabs:=2.19]
      #      params.by.cas.spec <- tmp[, .(CAS, Vdist,kelim, Rblood2plasma,
      #                               MW, hematocrit, million.cells.per.gliver)]
      #Since we don't yet have 2comp predictions,
      #just choose some arbitrary numbers as starting points for the fit.
      params.by.cas.spec[, V1:=1]
      params.by.cas.spec[, Ralphatokelim:=1.2]
      params.by.cas.spec[, Fbetaofalpha:=0.8]

      #      if (modelfun=="analytic") params.by.cas.spec[, setdiff(names(params.by.cas.spec),
      #                                                        c("CAS",
      #                                                          "kelim",
      #                                                          "Ralphatokelim",
      #                                                          "Fbetaofalpha",
      #                                                          "V1",
      #                                                          "Fgutabs",
      #                                                          "kgutabs")):=NULL]
    }

    ###########################################################################################################
    ###########################################################################################################
    ###########################################################################################################

    ### merge data.set and params.by.cas.spec to add parameter columns
    data.set <- merge(data.set,params.by.cas.spec,by=c("Species","Compound","CAS"))
    ### create vector of parameter names using column names of params.by.cas.spec
    paramnames <- names(params.by.cas.spec)
    ### remove names that were only used for merging, i.e. CAS, Species, and Compound
    paramnames <- paramnames[!(paramnames%in%c("CAS","Species","Compound"))]

    #Replace spaces in references with "."
    data.set[, Reference:=gsub(Reference,
                               pattern=' ',
                               replacement='.')]

    ### PK.fit.joint is a data.frame containing a row of parameter values per param.value.type per CAS
    PK.fit.joint <- data.set[,
                             analyze.pk.data(fitdata=.SD, ### what is .SD
                                             this.cas=CAS,
                                             paramnames=paramnames,
                                             modelfun=modelfun,
                                             model=model),
                             by=c("CAS","Species")]

    ### that correspond to multiple references
    multi.ref.cas <- data.set[,length(unique(Reference))>1,by=c("CAS","Species")]
    multi.ref.cas <- subset(multi.ref.cas,V1==T)$CAS
    if (length(multi.ref.cas)>0)
    {
      data.set.multi.ref <- subset(data.set,CAS %in% multi.ref.cas)


      PK.fit.separate <- data.set.multi.ref[,
                                            analyze.pk.data(fitdata=.SD,
                                                            this.cas=CAS,
                                                            paramnames=paramnames,
                                                            modelfun=modelfun,
                                                            model=model,
                                                            this.reference=Reference),
                                            by=c("CAS","Species","Reference")]



      PK.fit.bind <- rbind(PK.fit.joint,
                           PK.fit.separate)
    } else PK.fit.bind <- PK.fit.joint

    #record which model was fit to data and whether it was full or analytical
    PK.fit.bind[, model:=model]
    PK.fit.bind[, model.type:=modelfun]

    ###########################################################################################################
    ###########################################################################################################
    ###########################################################################################################

    #Do post fit calculations:
    ### basically calculate parameters based on previously calculated paramaters
    PK.fit.table <- PK.fit.bind
    if (model=="1compartment")
    {
      #Get stats for fitted total clearance :    L/kg body weight/h
      PK.fit.table[,CLtot:=Vdist*kelim]

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
                   Css:=ifelse(is.na(Fgutabs),1,Fgutabs)/(24*CLtot)] # 1 mg/kg/day / L/day/kg -> mg/L

      #Get statistics for halflife from fitted values
      PK.fit.table[,
                   halflife:=log(2)/kelim]


      PK.fit.table[,
                   tpeak.oral:=log(kgutabs/kelim)/(kgutabs-kelim)]

      ### why does Fgutabs have to be 1 right here
      ### this just copies the Ccompartment column
      ### not sure how values correspond to specific param.value.types
      ### again, it's just a time, conc, and auc matrix
      PK.fit.table[,
                   Cpeak.oral.1mgkg:=analytic_1comp_fun(
                     params=list(
                       Fgutabs = 1,
                       kgutabs = kgutabs,
                       kelim =kelim,
                       Vdist = Vdist
                     ),
                     dose=1,
                     tpeak.oral,
                     iv.dose=F)[,"Ccompartment"]]

      ### subset only certain param.value.types, not sure why
      PK.fit.table <- PK.fit.table[param.value.type%in%c("Predicted","Fitted geometric mean","Fitted geometric std dev")]

    } else if(model=="2compartment")
    {
      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   c("beta",
                     "alpha"):=list(Fbetaofalpha*Ralphatokelim*kelim,
                                    Ralphatokelim*(kelim+10^-6))]

      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   c("k21",
                     "k12"):=list(alpha*beta/kelim,
                                  alpha + beta - kelim - alpha*beta/kelim)]

      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   c("halflife",
                     "Vss",
                     "CLtot",
                     "Varea.or.Vbeta"):=list(log(2)/beta,
                                             V1*(k21+k12)/k21,
                                             V1*(k21+k12)/k21*beta,
                                             V1*kelim/beta)]
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
      PK.fit.table[param.value.type!="Predicted",
                   Css:=(Fgutabs*1)/(CLtot*24)]

      PK.fit.table <- PK.fit.table[param.value.type%in%c("Predicted","Fitted geometric mean","Fitted geometric std dev")]
    }
  }

  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################

  PK.fit.table <- PK.fit.table[order(Compound,Species)]

  ### figure out way without converting to and back from data.frame
  ### coerce data.set to data.frame
  data.set <- as.data.frame(data.set)

  ### function to apply loq fix to a single data.frame
  loq_fix <- function(data) {
    if(all(is.na(data$LOQ))) {data$LOQ <- 0.45 * min(data$Value, na.rm = TRUE)}
    return(data)
  }

  ### split data.set into list of data.frame objects
  split_df <- split(data.set, list(data.set$Compound, data.set$Reference, data.set$Media), drop = TRUE)
  # split_df <- split(data.set, data.set$Compound)

  ### apply loq_fix to each data.set
  split_df_loq <- lapply(split_df, loq_fix)

  ### 'unsplit' data.sets
  data.set <- do.call(rbind, split_df_loq)

  rownames(data.set) <- c()

  ### coerce data.set back to data.table
  data.set <- as.data.table(data.set)

  out <- list(PK.fit.table, data.set)
  return(out)
}



