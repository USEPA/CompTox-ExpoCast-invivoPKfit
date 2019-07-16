#' plot_conctime
#'
#' Creates a PDF of the fits from \code{\link{fit_all}} and data together
#'
#' @param PK.fit.table Table of parameter estimates from \code(\link(fit_all))
#' @param data.set A data.table of concentration vs. time data for a given
#'  chemical
#' @param model Analytic model to evaluate. Currently only "1compartment" or
#'  "2compartment" are implemented.
#' @param mean.type (Defaults to fitted geometric mean)
#' @param fit.factors Add a trend-line using the fit for each combination of the
#'  factors listed in this character vector. (defaults to Compound, 
#'  Dose.nominal, Data.Analyzed, and Route)
#' @param plot.split.factor For each compound, make a separate plot for each
#'  unique value of this factor (defaults to Route)
#' #param shape.factor Indicate this factor by shape in the plots (defaulta to
#'  Data.Analyzed)
#' @param color.factor Indicate this factor by color in the plots (defaults to
#'  Data.Analyzed')
#' @param Add a curve for the HTTK predictions, if available (default = FALSE)
#' @param omit.zero Drop zero dose values from plots (default = TRUE)
#'
#' @return None
#'
#' @import RColorBrewer
#' @importfrom scales scientific_format
#' @import ggplot2
#'
#' @export
#'
#' @author{Caroline Ring, John Wambaugh}
plot_conctime <- function(PK.fit.table,
                          data.set,
                          model,
                          mean.type="Fitted geometric mean",
                          fit.factors=c("Compound", 
                                    "Data.Analyzed", 
                                    "Dose.nominal", 
                                    "Route",
                                    "LOQ"),
                          plot.split.factor="Route",
                          shape.factor="Data.Analyzed",
                          color.factor="Data.Analyzed",
                          plot.httk.pred=F,
                          omit.zero=T,
                          plot.loq=T,
                          fit.plot.points=25)
{
  scientific_10 <- function(x) {                                  
    out <- gsub("1e", "10^", scales::scientific_format()(x))              
    out <- gsub("\\+","",out)                                     
    out <- gsub("10\\^01","10",out)                               
    out <- parse(text=gsub("10\\^00","1",out))                    
  }  

  data.set <- copy(data.set)
  if (omit.zero) data.set <- data.set[Dose>0,]
  
  if (plot.httk.pred)
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
           by=.(Compound, CAS, Source, Route, Dose)]
  data.set[, Dose.nominal.units.type:=paste(Dose.nominal, "mg/kg", Route)]

  setkey(data.set, Compound, CAS, Source)
  #Replace spaces in data.set references with periods, to match PK.fit.table style
  data.set[, Reference:=gsub(x=Reference, 
                             pattern=' ', 
                             replacement='.')]

  #create a nominal dose + route column for labeling plots later
  data.set[, studyid.nominal:=paste(Dose.nominal, 'mg/kg', Route)]
  
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
  if (model=="1compartment")
  {
    params.in <- list(Vdist=Vdist,
      kelim=kelim,
      kgutabs=kgutabs,
      Fgutabs=Fgutabs)
  } else if (model == "2compartment") {
    params.in<-list(Fbetaofalpha=Fbetaofalpha,
      Ralphatokelim=Ralphatokelim,
      kelim=kelim,
      V1=V1,
      kgutabs=kgutabs,
      Fgutabs=Fgutabs)
  } else stop()
    
# Select optimal colors
  colorvals <- RColorBrewer::brewer.pal(n=max(data.set[,
    length(unique(color.factor)),by=Compound]$V1),name='Set2')

  if (plot.httk.pred)
  {
    colorvals <- c('grey',colorvals)
    names(colorvals)[1] <- 'in vitro predicted'
  }
  
  pdf(paste(gsub(" ","",model),"-",gsub(" ","",mean.type),"-",Sys.Date(),".pdf",
    sep=""))
  for (cm in unique(data.set[, Compound])){ #plot each compound separately
  {
    this.compound.data <- data.set[Compound==cm]
    for (this.graph in sort(unique(this.compound.data[,
      plots.split.factor,with=F])))
    {
      this.graph.data <- this.compound.data[unlist(this.compound.data[,
        eval(plots.split.factor),with=F])==this.graph]
      this.graph.fits <- PK.fit.table[Compound=cm]
      this.graph.fits <- this.graph.fits[unlist(this.graph.fits[,
        eval(plots.split.factor),with=F])==this.graph]
      
      merge.table <- merge(this.graph.data,this.graph.fits,by=fit.factors)
      
# Determine x-axis range:
      min.time <- 0
      max.time <- round(max(unique(as.numeric(this.graph.data$Time)), 
        na.rm=T)*1.25)
        
# Pick points for evaluating the fits:
      plot.times <- seq(min.time,max.time,(max.time-min.time)/fit.plot.points)  

# Build a table of fit plots:
      bestfit.invivo<-merge.table[,evalfun(model=model,
             time=plot.times, 
             dose=unique(Dose.nominal), 
             route=unique(Route), 
             params.in = params.in),
             by=fit.factors]
      bestfit.invivo[, studyid.nominal:=paste(Dose.nominal, 'mg/kg', Route)]
      bestfit.invivo[Ccompartment<0,Ccompartment:=0]
      
# Determine y-axis range:
      min.conc <- 10^floor(log10(min(unique(as.numeric(this.graph.data$Value)),
        na.rm=T)))
      if (is.finite(mean(as.numeric(this.graph.data$LOQ)) &
        !is.na(mean(as.numeric(this.graph.data$LOQ))))) 
        min.conc <- min(min.conc,
        mean(as.numeric(this.graph.data))/5)
      max.conc <- 10^ceiling(log10(max(unique(as.numeric(this.graph.data$Value)),
        na.rm=T)))
  
      plotted.data <- sort(unique(this.graph.data[,(color.factor)]))
      if (plot.httk.pred)
      {
        names(colorvals)[2:(1+length(plotted.data))] <- plotted.data
      } else {
        names(colorvals) <- plotted.data
      }
      
      if (any(is.finite(this.graph.fits$AIC)))
      {
        p <- ggplot2::ggplot(data=data.set[Compound==cm, ]) +
          ggtitle(paste(cm," (",model,")",sep=""))+
          geom_line(data=bestfit.invivo,aes(x=time, 
            y=Ccompartment, 
            color=color.factor)) +
          geom_point(data=PK.fit.merge[Compound==cm & 
            Data.Analyzed!="Joint Analysis",], 
            aes(x=Time, y=Value, 
            color=color.factor, 
            shape=shape.factor), 
            size=3) +
          facet_wrap(facets=plots.split.factor, scales="free") +
          scale_x_continuous(limits = c(min.time,max.time)) +
          scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) + 
          scale_color_manual(values=colorvals) +
          xlab("Time (h)") + ylab("Concentration (mg/L)") + 
          theme_bw() +
          theme(aspect.ratio = 1)
        if (plot.loq)
        {
          p <- p + geom_hline(aes(yintercept=2*mean(LOQ)),color="Blue",
            linetype="dashed") +
            annotate("text",x=(max.time+min.time)/2,
            y=2*mean(bestfit.invivo[Compound==cm,"LOQ"]),
            label="Limit of Quantitiation")
        }
      } else {
        p <- ggplot2::ggplot(data=data.set[Compound==cm,]) +
          ggtitle(paste(cm," (",model,"): Optimizer Failed, No Curve Fit",sep=""))+
          geom_point(aes(x=Time, y=Value, color=color.factor, shape=shape.factor), size=3) +
          facet_wrap(facets=~studyid.nominal, scales="free") +
          geom_hline(aes(yintercept=2*mean(LOQ)),color="Blue",linetype="dashed")+
          annotate("text",x=(max.time+min.time)/2,y=2*mean(bestfit.invivo[Compound==cm,"LOQ"]),label="Limit of Quantitiation")+
          scale_x_continuous(limits = c(min.time,max.time)) +
          scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) + 
          scale_color_manual(values=colorvals) +
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