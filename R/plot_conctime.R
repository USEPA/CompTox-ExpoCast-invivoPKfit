#written by Caroline Ring
#modified by John Wambaugh

plot_conctime <- function(PK.fit.table,
                          data.set,
                          model,
                          mean.type="Fitted geometric mean")
{
  scientific_10 <- function(x) {                                  
    out <- gsub("1e", "10^", scientific_format()(x))              
    out <- gsub("\\+","",out)                                     
    out <- gsub("10\\^01","10",out)                               
    out <- parse(text=gsub("10\\^00","1",out))                    
  }  

  data.set <- copy(data.set)
  if (model=="1compartment")
  {
     PK.fit.table <- copy(PK.fit.table[param.value.type %in% c("Predicted",
                                                            mean.type), ])
  } else {
     PK.fit.table <- copy(PK.fit.table[param.value.type %in% mean.type, ])
  }
  
  
  #Get the nominal doses
  #get the nominal dose by rounding to 2 significant figure
  #(e.g. 5.11=5.1, 5.26=5.2, 1.01=1.0,0.98=1.0, etc.)
  data.set[, 
           Dose.nominal:=round(mean(Dose), digits=1),
           by=.(Compound, CAS, Source, Route)]
  data.set[, Dose.nominal.units.type:=paste(Dose.nominal, "mg/kg", Route)]
  setkey(data.set, Compound, CAS, Source)
  #Replace spaces in data.set references with periods, to match PK.fit.table style
  data.set[, Reference:=gsub(x=Reference, 
                             pattern=' ', 
                             replacement='.')]
  #Need to handle data from multiple labs
#  #So make a copy of the data for these two chemicals
#  data.joint <- data.set[Compound=='Bensulide' |
#                           Compound=='Propyzamide',]
#  #Change the Source for this data to Joint NHEERL/RTI
#  data.joint[, Source:='Joint NHEERL/RTI']
#  #Change the Reference to reflect the joint data
#  data.joint[, Reference:='RTI.2015, NHEERL.2015']
#  #Re-add this copy of the joint data to the data set.
#  data.set <- rbind(data.set, data.joint)
#  data.set[, Source:= factor(Source)]
  #merge together the data set and the in vivo PK fit results into one big data table
  data.set.joint <- subset(data.set,CAS %in% PK.fit.table[Data.Analyzed=="Joint Analysis",]$CAS)
  data.set.joint[,Reference:=paste(sort(unique(Reference)),collapse=", "),by=CAS]
  
  PK.fit.merge <- merge(rbind(data.set,data.set.joint),
                        PK.fit.table,
                        by=intersect(names(data.set), names(PK.fit.table)),
                        allow.cartesian=TRUE,all.y=T)
  
  #create a nominal dose + route column for labeling plots later
  PK.fit.merge[, studyid.nominal:=paste(Dose.nominal, 'mg/kg', Route)]
 # PK.fit.merge<-PK.fit.merge[!is.na(Time),]
  
  #Set up function for computing the best-fit line:
  #evaluates 1-compartment model using the params provided as arguments,
  #at the vector of times provided as the "time" argument,
  #with the dose and route provided as arguments
  evalfun <- function(model, time, dose, route,
                     params.in){
    
    params.in <- lapply(params.in, unique)
    
    out <- invivoPKfit:::analytic_model_fun(params=params.in, 
                                            dose=dose, 
                                            times=time, 
                                            time.units='h', 
                                            iv.dose=ifelse(route=='iv', TRUE, FALSE),
                                            model=model)
    return(as.data.table(out[, c('time', 'Ccompartment')]))
  }
  
  #get best-fit lines from in vivo fit parameter values
  if (model=="1compartment"){
  bestfit.invivo<-PK.fit.merge[, 
                               evalfun(model=unique(model),
                                       time=seq(from=0,
                                                to=max(Time,na.rm=T),
                                                length.out=100), 
                                       dose=unique(Dose.nominal), 
                                       route=unique(Route), 
                                       params.in = list(Vdist=Vdist,
                                                        kelim=kelim,
                                                        kgutabs=kgutabs,
                                                        Fgutabs=Fgutabs)), 
                               by=.(Compound, #eval model separately for each combination of chemical, source lab, nominal dose, and route
                                    Data.Analyzed, 
                                    Dose.nominal, 
                                    Route,
                                    param.value.type,
                                    LOQ)]
  } else if (model == "2compartment"){
    bestfit.invivo<-PK.fit.merge[, 
                                 evalfun(model=unique(model),
                                         time=seq(from=0,
                                                  to=max(Time,na.rm=T),
                                                  length.out=100), 
                                         dose=unique(Dose.nominal), 
                                         route=unique(Route), 
                                        params.in=list(Fbetaofalpha=Fbetaofalpha,
                                                    Ralphatokelim=Ralphatokelim,
                                                    kelim=kelim,
                                                    V1=V1,
                                                    kgutabs=kgutabs,
                                                    Fgutabs=Fgutabs)), 
                                 by=.(Compound, #eval model separately for each combination of chemical, source lab, nominal dose, and route
                                      Data.Analyzed, 
                                      Dose.nominal, 
                                      Route,
                                      param.value.type,
                                      LOQ)]
  }
  bestfit.invivo[, studyid.nominal:=paste(Dose.nominal, 'mg/kg', Route)]
  bestfit.invivo[param.value.type=="Predicted", Data.Analyzed:="in vitro predicted"]
  bestfit.invivo[Ccompartment<0,Ccompartment:=0]
  
  num.cols <-max(PK.fit.merge[,length(unique(Data.Analyzed)),by=Compound]$V1)
  colvals <- brewer.pal(n=num.cols,name='Set2')
  shapevals <- rep(21,length(colvals))
  if (model=="1compartment")
  {
    colvals <- c('grey',colvals)
    names(colvals)[1] <- 'in vitro predicted'
    shapevals <- c(NA, shapevals)
  }
  
  pdf(paste(gsub(" ","",model),"-",gsub(" ","",mean.type),"-",Sys.Date(),".pdf",sep=""))
  for (cm in unique(PK.fit.merge[, Compound])){ #plot each compound separately
  {
    min.time <- 0
    max.time <- round(max(unique(as.numeric(subset(PK.fit.merge[Compound==cm,],!is.na(Value))$Time)),na.rm=T)*1.25)
    min.conc <- 10^floor(log10(min(unique(as.numeric(PK.fit.merge[Compound==cm,]$Value)),na.rm=T)))
    
    if (is.finite(mean(as.numeric(PK.fit.merge[Compound==cm,]$LOQ))&!is.na(mean(as.numeric(PK.fit.merge[Compound==cm,]$LOQ))))) min.conc <- min(min.conc,mean(as.numeric(PK.fit.merge[Compound==cm,]$LOQ))/5)
    
    max.conc <- 10^ceiling(log10(max(unique(as.numeric(PK.fit.merge[Compound==cm,]$Value)),na.rm=T)))

    num.cols <- length(unique(PK.fit.merge[Compound==cm,Data.Analyzed]))
    plotted.data <- sort(unique(PK.fit.merge[Compound==cm,Data.Analyzed]))
    if (model=="1compartment")
    {
      
      names(colvals)[2:(1+length(plotted.data))] <- plotted.data
    } else {
      names(colvals) <- plotted.data
    }
    names(shapevals) <- names(colvals)
 #   if ("Joint Analysis" %in% names(shapevals)) shapevals["Joint Analysis"] <- NA
    
    if (is.finite(PK.fit.merge[Compound==cm,]$AIC[1]))
    {
      
      p <- ggplot(data=bestfit.invivo[Compound==cm, ]) +
        ggtitle(paste(cm," (",model,")",sep=""))+
        geom_line(aes(x=time, y=Ccompartment, color=Data.Analyzed)) +
        geom_point(data=PK.fit.merge[Compound==cm&Data.Analyzed!="Joint Analysis",], 
                   aes(x=Time, y=Value, color=Data.Analyzed, shape=Data.Analyzed), size=3) +
        facet_wrap(facets=~studyid.nominal, scales="free") +
        geom_hline(aes(yintercept=2*mean(LOQ)),color="Blue",linetype="dashed")+
        annotate("text",x=(max.time+min.time)/2,y=2*mean(bestfit.invivo[Compound==cm,"LOQ"]),label="Limit of Quantitiation")+
        scale_x_continuous(limits = c(min.time,max.time)) +
        scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) + 
        scale_color_manual(values=colvals) +
        scale_shape_manual(values=shapevals,guide=F)+
        xlab("Time (h)") + ylab("Concentration (mg/L)") + 
        theme_bw() +
        theme(aspect.ratio = 1)
    } else {
      p <- ggplot(data=PK.fit.merge[Compound==cm,]) +
        ggtitle(paste(cm," (",model,"): Optimizer Failed, No Curve Fit",sep=""))+
        geom_point(aes(x=Time, y=Value, color=Data.Analyzed, shape=Data.Analyzed), size=3) +
        facet_wrap(facets=~studyid.nominal, scales="free") +
        geom_hline(aes(yintercept=2*mean(LOQ)),color="Blue",linetype="dashed")+
        annotate("text",x=(max.time+min.time)/2,y=2*mean(bestfit.invivo[Compound==cm,"LOQ"]),label="Limit of Quantitiation")+
        scale_x_continuous(limits = c(min.time,max.time)) +
        scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) + 
        scale_color_manual(values=colvals) +
        scale_shape_manual(values=shapevals,guide=F)+
        xlab("Time (h)") + ylab("Concentration (mg/L)") + 
        theme_bw() +
        theme(aspect.ratio = 1) 
    }
  }
    
      print(p)
  }
  dev.off()
  return(0)
}