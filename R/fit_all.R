#written by Caroline Ring
#modified by John Wambaugh
#
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
#'
#' @return A data.table of fitted parameter values for each chemical.
#'
#' @export

fit_all <- function(data.set,
                    model,
                    modelfun=NA)
{

  data.set <- data.table::copy(data.set)

  
  #Ignore data close to LOQ:
  data.set[Value<2*LOQ,Value:=NA]

  #set an iv variable to TRUE/FALSE
  data.set[Route=='iv', iv:=TRUE]
  data.set[Route!='iv', iv:=FALSE]
  data.set[, Time.Days:=Time/24]  #convert time from hours to days
  data.set[, c('Max.Time.Days',
               'Time.Steps.PerHour'):=list(max(Time.Days),
                                        1/min(diff(c(0,sort(unique(Time)))))),
           by=.(CAS, Dose, Route)]
  
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
    params.by.cas <- data.table(CAS=sort(unique(data.set$CAS)),stringsAsFactors=F)
    #get the rat parameters for the 1-compartment model for each chemical
    if (model=='1compartment')
    {
      if (any(params.by.cas$CAS %in% get_cheminfo()))
        params.by.cas <- params.by.cas[CAS %in% get_cheminfo(), 
        httk::parameterize_1comp(chem.cas=CAS,
        default.to.human=TRUE,
        species="Rat"),
        by=CAS]
# Wambaugh et al. (2018) medians:
# apply(chem.invivo.PK.aggregate.data,2,function(x) median(as.numeric(x),na.rm=T))        
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), kelim:=0.25]
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), Vdist:=5.56]
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), Fgutabs:=0.31]
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), kgutabs:=2.19]
      
      if (modelfun=="analytic") params.by.cas[, setdiff(names(params.by.cas),
                                                        c("CAS",
                                                          "kelim",
                                                          "Vdist",
                                                          "Fgutabs",
                                                          "kgutabs")):=NULL]
      params.by.cas[, Fgutabs:=0.99] #Assume 99% bioavailability for everything
      params.by.cas[, kgutabs:=1]
    } else if (model=='2compartment') {
    # Use this when parameterize_2comp is implemented in httk
    #   params.by.cas <- data.set[, httk::parameterize_2comp(chem.cas=CAS,
    #                                                  default.to.human=TRUE,
    #                                                  species="Rat"),
    #                             by=CAS]
      if (any(params.by.cas$CAS %in% get_cheminfo()))
        params.by.cas <- params.by.cas[CAS %in% get_cheminfo(), 
        httk::parameterize_1comp(chem.cas=CAS,
        default.to.human=TRUE,
        species="Rat"),
        by=CAS]
# Wambaugh et al. (2018) medians:
# apply(chem.invivo.PK.aggregate.data,2,function(x) median(as.numeric(x),na.rm=T))        
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), kelim:=0.25]
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), Vdist:=5.56]
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), Fgutabs:=0.31]
      params.by.cas[!params.by.cas$CAS%in%get_cheminfo(), kgutabs:=2.19]
      params.by.cas <- tmp[, .(CAS, Vdist,kelim, Rblood2plasma,
                               MW, hematocrit, million.cells.per.gliver)]
      #Since we don't yet have 2comp predictions,
      #just choose some arbitrary numbers as starting points for the fit.
      params.by.cas[, V1:=1]
      params.by.cas[, Ralphatokelim:=1.2]
      params.by.cas[, Fbetaofalpha:=0.8]
  #    params.by.cas[, kelim:=max(-lm(log(Value)~Time,subset(chem.invivo.PK.data,Compound=="Carbendazim"))$coefficients[2],10^-5)]
      params.by.cas[, Fgutabs:=0.99] #Assume 99% bioavailability for everything
      params.by.cas[, kgutabs:=1]
    
      if (modelfun=="analytic") params.by.cas[, setdiff(names(params.by.cas),
                                                        c("CAS",
                                                          "kelim",
                                                          "Ralphatokelim",
                                                          "Fbetaofalpha",
                                                          "V1",
                                                          "Fgutabs",
                                                          "kgutabs")):=NULL]
    }
    data.set <- merge(data.set,params.by.cas,by="CAS")
    paramnames <- names(params.by.cas)
    paramnames <- paramnames[paramnames!='CAS']
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
                             by=CAS]

    multi.ref.cas <- data.set[,length(unique(Reference))>1,by=CAS]
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
                               by=c("CAS","Reference")]
                               
  
      
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
      
      #Get HTTK-predicted total clearance for each CAS    L/h/kg BW
      PK.fit.table[param.value.type=="Predicted",
                   CLtot:=httk::calc_total_clearance(chem.cas=CAS,
                                               species="Rat",
                                               default.to.human=TRUE,
                                               suppress.messages=TRUE),
                   by=CAS]
  
      #Get stats for fitted total clearance :    L/kg body weight/h
      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   CLtot:=Vdist*kelim]
  #    PK.fit.table[param.value.type %in% c("Fitted arithmetic std dev",
  #                                         "Fitted geometric std dev"),
  #                 CLtot:=sqrt(Vdist^2+kelim^2)]
  
      #Get HTTK-predicted Css for 1-compartment model
      PK.fit.table[param.value.type=="Predicted",
                   Css:=httk::calc_analytic_css(chem.cas=CAS,
                                          species="Rat",
                                          output.units='mg/L',
                                          model='1compartment',
                                          default.to.human=TRUE,
                                          suppress.messages=TRUE),
                   by=CAS]
  
      #Get statistics for Css from fitted CLtot values
      #Get stats for fitted total clearance:
      PK.fit.table[param.value.type %in% c("Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   Css:=Fgutabs/(24*CLtot)] # 1 mg/kg/day / L/day/kg -> mg/L
  
      # PK.fit.table[, Css.arith.se := -Css.arith* #propagation of uncertainty
      #            CLtot.arith.se/
      #            CLtot.arith]
      # PK.fit.table[, Css.geo.se := -Css.geo* #propagation of uncertainty
      #            CLtot.geo.se/
      #            CLtot.geo]
  
      #Get statistics for halflife from fitted values
      PK.fit.table[param.value.type %in% c("Predicted",
                                           "Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   halflife:=log(2)/kelim]
  
  
      PK.fit.table[param.value.type %in% c("Predicted",
                                           "Fitted arithmetic mean",
                                           "Fitted geometric mean",
                                           "Fitted mode"),
                   tpeak.oral:=log(kgutabs/kelim)/(kgutabs-kelim)]
     
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
      PK.fit.table[param.value.type=="Predicted",
                   Css:=httk::calc_analytic_css(chem.cas=CAS,
                                                species="Rat",
                                                output.units='mg/L',
                                                model='1compartment',
                                                default.to.human=TRUE,
                                                suppress.messages=TRUE),
                   by=CAS]
    }
  }
  
  return(PK.fit.table)
}


