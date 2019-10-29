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
#' #param ratio between the weights used to report the data and the weights used 
#'   for the dose. For example, ug/L data and mg/kg/day dose would be 0.001
#'   (defaults to 1) 
#'
#' @return A data.table of fitted parameter values for each chemical.
#'
#' @author Caroline Ring, John Wambaugh
#'
#' @export

fit_all <- function(data.set,
                    model,
                    modelfun=NA,
                    ratio.data.to.dose=1)
{

  data.set <- data.table::copy(data.set)
  
  N.PREV <- dim(data.set)[1]
  cat(paste(N.PREV,"concentration vs. time observations loaded.\n"))
  cat(paste(length(unique(data.set$CAS)),"unique chemicals,",
        length(unique(data.set$Species)),"unique species, and",
        length(unique(data.set$Reference)),"unique references.\n"))
  
  if (is.character(class(data.set$Dose)))
  {
    cat("Column \"Dose\" converted to numeric.\n")
    data.set$Dose <- as.numeric(data.set$Dose)
  }
  # Right now code only recognizes "po" and "iv" as routes (my bad):
  data.set[Route=="oral",Route:="po"]
  data.set[Route=="intravenous",Route:="iv"]
  data.set <- data.set[Route %in% c("po","iv")]
  cat(paste("Restricting to intravenous and oral routes eliminates",
      N.PREV - dim(data.set)[1],"observations.\n"))
  cat(paste(dim(data.set)[1],"observations of",
        length(unique(data.set$CAS)),"unique chemicals,",
        length(unique(data.set$Species)),"unique species, and",
        length(unique(data.set$Reference)),"unique references remain.\n"))
  N.PREV <- dim(data.set)[1]
  
  # Harmonize the compound names:
  data.set[,Compound:=tolower(Compound)]
  
  # This way the weight units cancel (must still pay attention to denominator
  # of data to determine units for Vd):
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
  
  # Because doses are administered instantaneously in the model, we can't 
  # handle early time points below LOQ:
  data.set[,
           Usable:=ifelse(Time >= .SD[Value==max(Value,na.rm=T),Time] |
                 !is.na(Value),T,F),
           by=.(Route,Reference,CAS,Species)]
  data.set <- data.set[Usable==TRUE]
  cat(paste("Eliminating observations for doses that are below LOQ before the peak conc. is reached eliminates",
      N.PREV - dim(data.set)[1],"observations.\n"))
  cat(paste(dim(data.set)[1],"observations of",
        length(unique(data.set$CAS)),"unique chemicals,",
        length(unique(data.set$Species)),"unique species, and",
        length(unique(data.set$Reference)),"unique references remain.\n"))
  N.PREV <- dim(data.set)[1]

  #Non-comapartmental fits:
  if (model=="noncompartment")
  {
    # Normalize all concentrations by dose to fit jointly:
    data.set[,Value.Norm:=Value/Dose]

    PK.fit.table<- NULL
    for (this.cas in sort(unique(data.set$CAS)))
    {
      this.subset <- subset(data.set,CAS==this.cas&!is.na(Value.Norm))
             
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
                      Mean.iv.Dose=mean(subset(this.subset,Route=="iv")$Dose),
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
          if (any(is.na(this.route.subset$Subject)))
          {
            for (this.reference in unique(this.route.subset$Reference))
            {
              for (this.dose in unique(this.route.subset$Dose))
              {
                this.route.subset[is.na(this.route.subset$Subject)&this.route.subset$Dose==this.dose,"Subject"] <- paste(this.reference,this.route,this.dose,sep="-")
              }
            }
          }
          this.route.subset$id <- paste(this.row[1,"Reference"],this.route.subset$Subject,sep="-")
          this.route.subset$conc <- this.route.subset$Value.Norm
          this.route.subset$time <- this.route.subset$Time
          for (this.id in unique(this.route.subset$id))
          {
            if (dim(subset(this.route.subset,id==this.id))[1]<3)
            {
              this.route.subset <- subset(this.route.subset,id!=this.id)
            }
          }
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
            if (length(this.route.subset$time)==length(unique(this.route.subset$time)))
            {
              out <- try(nca.complete(data=this.route.subset,dose=dose.arg,method="z"))
            }
            else
            {
              out <- try(nca.batch(data=this.route.subset,dose=dose.arg,method="z"))
            }
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
        PK.fit.table <- rbind(PK.fit.table,this.row)
      }
    }
    
    PK.fit.table <- as.data.table(PK.fit.table)
    
    
    # Multiply by dose to return from normalized uinits:
    PK.fit.table[,Vd.iv:=Vd.iv*Mean.iv.Dose]
    PK.fit.table[,Vd.po:=Vd.po*Mean.po.Dose]

  } else {
  
    params.by.cas.spec <- data.set[,unique(.SD[,.(Compound)]),by=.(CAS,Species)]

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
        params.by.cas.spec[CAS %in% get_cheminfo(),
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
    data.set <- merge(data.set,params.by.cas.spec,by=c("Species","Compound","CAS"))
    paramnames <- names(params.by.cas.spec)
    paramnames <- paramnames[!(paramnames%in%c("CAS","Species","Compound"))]
    #Replace spaces in references with "."
    data.set[, Reference:=gsub(Reference,
                            pattern=' ',
                            replacement='.')]

    PK.fit.joint <- data.set[,             
                             analyze.pk.data(fitdata=.SD,
                                             this.cas=CAS,
                                             paramnames=paramnames,
                                             modelfun=modelfun,
                                             model=model),
                             by=c("CAS","Species")]

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
    
                  
    #Do post fit calculations:
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
      PK.fit.table[,
                   Css:=ifelse(is.na(Fgutabs),1,Fgutabs)/(24*CLtot)] # 1 mg/kg/day / L/day/kg -> mg/L
    
      #Get statistics for halflife from fitted values
      PK.fit.table[,
                   halflife:=log(2)/kelim]
  
  
      PK.fit.table[,
                   tpeak.oral:=log(kgutabs/kelim)/(kgutabs-kelim)]
 
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

    }
    }
  
  PK.fit.table <- PK.fit.table[order(Compound,Species)]
  
  return(PK.fit.table[param.value.type%in%c("Predicted","Fitted geometric mean","Fitted geometric std dev")])
}


