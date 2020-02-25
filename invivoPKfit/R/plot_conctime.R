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
  #All the TNO Lit. Review chemicals are already listed as nominal doses
  data.set[Source=="TNO Lit. Review", Dose.nominal:=Dose]
  #For RTI and NHEERl chemicals:
  #get the nominal dose by rounding to 1 significant figure
  #(e.g. 5.11=5, 5.26=5, 1.01=1,0.98=1, etc.)
  data.set[Source!="TNO Lit. Review", 
           Dose.nominal:=round(mean(Dose), digits=1),
           by=.(Compound, CAS, Source, Route)]
  data.set[Source=="RTI 2015" & Route=="po", Dose.nominal:=round(Dose.nominal)]
  data.set[, Dose.nominal.units.type:=paste(Dose.nominal, "mg/kg", Route)]
  setkey(data.set, Compound, CAS, Source)
  #Replace spaces in data.set references with periods, to match PK.fit.table style
  data.set[, Reference:=gsub(x=Reference, 
                             pattern=' ', 
                             replacement='.')]
  #Need to handle joint RTI/NHEERL data
  #Know already that the two shared chemicals are bensulide and propyzamide
  #So make a copy of the data for these two chemicals
  data.joint <- data.set[Compound=='Bensulide' |
                           Compound=='Propyzamide',]
  #Change the Source for this data to Joint NHEERL/RTI
  data.joint[, Source:='Joint NHEERL/RTI']
  #Change the Reference to reflect the joint data
  data.joint[, Reference:='RTI.2015, NHEERL.2015']
  #Re-add this copy of the joint data to the data set.
  data.set <- rbind(data.set, data.joint)
  data.set[, Source:= factor(Source)]
  #merge together the data set and the in vivo PK fit results into one big data table
  PK.fit.merge <- merge(data.set,
                        PK.fit.table,
                        by=intersect(names(data.set), names(PK.fit.table)),
                        allow.cartesian=TRUE)
  #create a nominal dose + route column for labeling plots later
  PK.fit.merge[, studyid.nominal:=paste(Dose.nominal, 'mg/kg', Route)]
  
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
                                       time=seq(from=0, #eval every second, for smooth line
                                                to=max(Time),
                                                length.out=max(Time)*60), 
                                       dose=unique(Dose.nominal), 
                                       route=unique(Route), 
                                       params.in = list(Vdist=Vdist,
                                                        kelim=kelim,
                                                        kgutabs=kgutabs,
                                                        Fgutabs=Fgutabs)), 
                               by=.(Compound, #eval model separately for each combination of chemical, source lab, nominal dose, and route
                                    Source, 
                                    Dose.nominal, 
                                    Route,
                                    param.value.type,
                                    LOQ)]
  } else if (model == "2compartment"){
    bestfit.invivo<-PK.fit.merge[, 
                                 evalfun(model=unique(model),
                                         time=seq(from=0,
                                                  to=max(Time),
                                                  length.out=max(Time)*60), 
                                         dose=unique(Dose.nominal), 
                                         route=unique(Route), 
                                        params.in=list(Fbetaofalpha=Fbetaofalpha,
                                                    Ralphatokelim=Ralphatokelim,
                                                    kelim=kelim,
                                                    V1=V1,
                                                    kgutabs=kgutabs,
                                                    Fgutabs=Fgutabs)), 
                                 by=.(Compound, #eval model separately for each combination of chemical, source lab, nominal dose, and route
                                      Source, 
                                      Dose.nominal, 
                                      Route,
                                      param.value.type,
                                      LOQ)]
  }
  bestfit.invivo[, studyid.nominal:=paste(Dose.nominal, 'mg/kg', Route)]
  bestfit.invivo[param.value.type=="Predicted", Source:="in vitro predicted"]
  bestfit.invivo[Ccompartment<0,Ccompartment:=0]
  
  if (model=="1compartment")
  {
    #Set up vector of colors to denote source lab in plots
    colvals <- c('grey',brewer.pal(n=4,name='Set2'))
    names(colvals) <- c('in vitro predicted',
                        'Joint NHEERL/RTI',
                        'NHEERL 2015',
                        'RTI 2015',
                        'TNO Lit. Review')
    
    shapevals <- c(NA, NA, 21, 21, 21)
    names(shapevals) <- c('in vitro predicted',
                          'Joint NHEERL/RTI',
                          'NHEERL 2015',
                          'RTI 2015',
                          'TNO Lit. Review')
  } else {
    #Set up vector of colors to denote source lab in plots
    colvals <- brewer.pal(n=4,name='Set2')
    names(colvals) <- c('Joint NHEERL/RTI',
                        'NHEERL 2015',
                        'RTI 2015',
                        'TNO Lit. Review')
    
    shapevals <- c(NA, 21, 21, 21)
    names(shapevals) <- c('Joint NHEERL/RTI',
                          'NHEERL 2015',
                          'RTI 2015',
                          'TNO Lit. Review')
  }
  
  pdf(paste(gsub(" ","",model),"-",gsub(" ","",mean.type),"-",Sys.Date(),".pdf",sep=""))
  for (cm in unique(PK.fit.merge[, Compound])){ #plot each compound separately
  {
    min.time <- 0
    max.time <- round(max(unique(as.numeric(subset(PK.fit.merge[Compound==cm,],!is.na(Value))$Time)))*1.25)
    min.conc <- 10^floor(log10(min(unique(as.numeric(PK.fit.merge[Compound==cm,]$Value)),na.rm=T)))
    if (is.finite(mean(as.numeric(PK.fit.merge[Compound==cm,]$LOQ))&!is.na(mean(as.numeric(PK.fit.merge[Compound==cm,]$LOQ))))) min.conc <- min(min.conc,mean(as.numeric(PK.fit.merge[Compound==cm,]$LOQ))/5)
    max.conc <- 10^ceiling(log10(max(unique(as.numeric(PK.fit.merge[Compound==cm,]$Value)),na.rm=T)))
    if (is.finite(PK.fit.merge[Compound==cm,]$AIC[1]))
    {
      
      p <- ggplot(data=bestfit.invivo[Compound==cm, ]) +
        ggtitle(paste(cm," (",model,")",sep=""))+
        geom_line(aes(x=time, y=Ccompartment, color=Source)) +
        geom_point(data=PK.fit.merge[Compound==cm,], 
                   aes(x=Time, y=Value, color=Source, shape=Source), size=3) +
        facet_wrap(facets=~studyid.nominal, scales="free") +
        geom_hline(aes(yintercept=2*mean(LOQ)),color="Blue",linetype="dashed")+
        scale_x_continuous(limits = c(min.time,max.time)) +
        scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) + 
        scale_color_manual(values=colvals) +
        scale_shape_manual(values=c("Joint NHEERL/RTI"=NA,
                                    "NHEERL 2015"=21,
                                    "RTI 2015"=21,
                                    "TNO Lit. Review"=21))+
        xlab("Time (h)") + ylab("Concentration (mg/L)") + 
        theme_bw() +
        theme(aspect.ratio = 1)
    } else {
      p <- ggplot(data=PK.fit.merge[Compound==cm,]) +
        ggtitle(paste(cm," (",model,"): Optimizer Failed, No Curve Fit",sep=""))+
        geom_point(aes(x=Time, y=Value, color=Source, shape=Source), size=3) +
        facet_wrap(facets=~studyid.nominal, scales="free") +
        geom_hline(aes(yintercept=2*mean(LOQ)),color="Blue",linetype="dashed")+
        scale_x_continuous(limits = c(min.time,max.time)) +
        scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) + 
        scale_color_manual(values=colvals) +
        scale_shape_manual(values=c("Joint NHEERL/RTI"=NA,
                                    "NHEERL 2015"=21,
                                    "RTI 2015"=21,
                                    "TNO Lit. Review"=21))+
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