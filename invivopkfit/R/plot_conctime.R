#' plot_conctime
#'
#' Creates a PDF of the fits from \code{\link{fit_all}} and data together
#'
#'
#' @param PK.fit.table Table of parameter estimates from \code{\link{fit_all}}
#' @param data.set A data.table of concentration vs. time data for a given
#' chemical
#' @param model Analytic model to evaluate. Currently only "1compartment" or
#' "2compartment" are implemented.
#' @param mean.type (Defaults to fitted geometric mean)
#' @param fit.factors Add a trend-line using the fit for each combination of
#' the factors listed in this character vector. (defaults to Compound,
#' Dose.nominal, Data.Analyzed, and Route)
#' @param shape.factor placeholder
#' @param plot.split.factor For each compound, make a separate plot for each
#' unique value of this factor (defaults to Route) #param shape.factor Indicate
#' this factor by shape in the plots (defaulta to Data.Analyzed)
#' @param color.factor Indicate this factor by color in the plots (defaults to
#' Data.Analyzed')
#' @param plot.httk.pred placeholder
#' @param omit.zero Drop zero dose values from plots (default = TRUE)
#' @param plot.loq placeholder
#' @param Add a curve for the HTTK predictions, if available (default = FALSE)
#' @param fit.plot.points placeholder
#' @return None
#' @author Caroline Ring, John Wambaugh
#' @export plot_conctime
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
                          # shape.factor="Data.Analyzed",
                          color.factor="Data.Analyzed",
                          plot.httk.pred=F,
                          omit.zero=T,
                          plot.loq=T,
                          fit.plot.points=200)
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
  data.set[, Dose.nominal.units.type:=paste(
    paste("Dose:", Dose.nominal, "mg/kg"),
    paste("Route:", Route),
    paste("Reference:", Reference),
    paste("Media:", Media),
    paste("LOQ:", LOQ), sep = "\n")]

  setkey(data.set, Compound, CAS, Source)
  #Replace spaces in data.set references with periods, to match PK.fit.table style
  data.set[, Reference:=gsub(x=Reference,
                             pattern=' ',
                             replacement='.')]

  #create a nominal dose + route column for labeling plots later
  data.set[, studyid.nominal:=paste(Dose.nominal, 'mg/kg', Route)]

  ### make all compound names completely lower-case
  data.set[,Compound:=tolower(Compound)]

  #Set up function for computing the best-fit line:
  #evaluates 1-compartment model using the params provided as arguments,
  #at the vector of times provided as the "time" argument,
  #with the dose and route provided as arguments
  evalfun <- function(model, time, dose, route, reference, media, loq, params.in) ###added reference, media, and loq
  {
    params.in <- lapply(params.in, unique)

    out <- invivoPKfit:::analytic_model_fun(params=params.in,
                                            dose=dose,
                                            times=time,
                                            time.units='h',
                                            iv.dose=ifelse(route=='iv', TRUE, FALSE),
                                            model=model)

    out <- as.data.table(out[, c('time', 'Ccompartment')])
    out[,Dose.nominal:=dose]
    out[,Route:=route]
    out[,Model:=model]
    out[,Media:=media] ### added
    out[,Reference:=reference] ### added
    out[,LOQ:=loq] ### added

    return(out)
  }

  # Select optimal colors
  colorvals <- RColorBrewer::brewer.pal(n=max(data.set[,
                                                       length(unique(color.factor)),by=Compound]$V1),name='Set2')

  if (plot.httk.pred)
  {
    colorvals <- c('grey',colorvals)
    names(colorvals)[1] <- 'in vitro predicted'
  }

  grDevices::pdf(paste(gsub(" ","",model),"-",gsub(" ","",mean.type),"-",Sys.Date(),".pdf",
                       sep=""))

  for (this.compound in unique(data.set[, Compound])) #plot each compound separately
  {
    this.compound.data <- data.set[Compound==this.compound]

    this.reference.list <- list()
    for (this.reference in unique(this.compound.data[, Reference])) {

      ### coerce NA "Value" to 0.5 * mean LOQ
      this.reference.data <- this.compound.data[Reference == this.reference]

      this.medium.list <- list()
      for(this.medium in unique(this.reference.data[, Media])) {

        this.medium.data <- this.reference.data[Media == this.medium]

        this.medium.data[, "Value"][is.na(this.medium.data[, "Value"])] <- 0.5 * this.medium.data[, mean(LOQ)]

        ### do this without converting to and back from a data.frame
        this.medium.data <- as.data.frame(this.medium.data)
        this.medium.data$L.O.Q. <- ifelse(this.medium.data$Value <= mean(this.medium.data$LOQ), "Below", "Above")
        this.medium.data <- as.data.table(this.medium.data)
        this.medium.data <- this.medium.data[, LOQ. := as.factor(L.O.Q.)]

        this.medium.list[[this.medium]] <- this.medium.data
      }
      this.reference.data <- do.call(rbind, this.medium.list)

      this.reference.list[[this.reference]] <- this.reference.data
    }

    this.compound.data <- do.call(rbind, this.reference.list)

    ### create new binary column to determine shape for plots
    ### where values less than or equal to mean loq are assigned '1'
    ### and values greater than mean loq are assigned '0'
    # this.compound.data <- as.data.frame(this.compound.data)
    # this.compound.data$shape.factor <- ifelse(this.compound.data$Value <= mean(this.compound.data$LOQ), "Below", "Above")
    # this.compound.data <- as.data.table(this.compound.data)
    # this.compoound.data <- this.compound.data[, shape.factor := as.factor(shape.factor)]

    print(this.compound)

    this.graph.fits <- PK.fit.table[Compound==this.compound]

    this.dose.regimens <- this.compound.data[,(unique(c("Dose","Route", "Media", "Reference", ### added media
                                                        intersect(fit.factors,colnames(this.compound.data))))),with=F]
    this.dose.regimens <- this.dose.regimens[!duplicated(this.dose.regimens)]

    this.dose.regimens[,Compound:=tolower(Compound)]

    merge.table <- merge(this.dose.regimens,this.graph.fits,
                         by=intersect(intersect(fit.factors,
                                                names(this.graph.fits)),
                                      names(this.compound.data)), allow.cartesian=T)

    merge.table[,Reference.y:=NULL] ### delete duplicate Reference column
    colnames(merge.table)[colnames(merge.table) == "Reference.x"] <- "Reference" ### Rename column

    # Determine x-axis range:
    min.time <- 0
    max.time <- round(max(unique(as.numeric(this.compound.data$Time)),
                          na.rm=T)*1.25)

    # Pick points for evaluating the fits:
    plot.times <- seq(min.time,max.time,(max.time-min.time)/fit.plot.points)

    ### if merge.table includes data from multiple references, only include that data for making plots
    ### otherwise, each reference specific model will be plotted
    if(any(merge.table$Data.Analyzed == "Joint Analysis")) {
      merge.table <- merge.table[Data.Analyzed == "Joint Analysis"]

      ### do this without converting to and back from data.frame
      merge.table <- as.data.frame(merge.table)

      merge.table <- separate_rows(merge.table, Reference, sep = ", ")

      merge.table <- semi_join(merge.table, data.set, by = c("Dose", "Route", "Media", "Reference"))

      merge.table <- as.data.table(merge.table)
    }

    # Build a table of fit plots:
    #get best-fit lines from in vivo fit parameter values
    if (model=="1compartment")
    {
      bestfit.invivo<-merge.table[, evalfun(model="1compartment",
                                            time=plot.times,
                                            dose=Dose.nominal,
                                            route=Route,
                                            media = Media, ### added
                                            reference = Reference, ### added
                                            loq = LOQ, ### added
                                            params.in = list(Vdist=Vdist,
                                                             kelim=kelim,
                                                             kgutabs=kgutabs,
                                                             Fgutabs=Fgutabs)),
                                  by = seq_len(nrow(merge.table))]
    } else if (model == "2compartment") {
      bestfit.invivo<-merge.table[, evalfun(model=model,
                                            time=plot.times,
                                            dose=Dose.nominal,
                                            route=Route,
                                            media = Media, ### added
                                            reference = Reference, ### added
                                            loq = LOQ, ### added
                                            params.in = list(Fbetaofalpha=Fbetaofalpha,
                                                             Ralphatokelim=Ralphatokelim,
                                                             kelim=kelim,
                                                             V1=V1,
                                                             kgutabs=kgutabs,
                                                             Fgutabs=Fgutabs)),
                                  by = seq_len(nrow(merge.table))]
    } else stop()
    bestfit.invivo[, Dose.nominal.units.type:=paste(
      paste("Dose:", Dose.nominal, "mg/kg"),
      paste("Route:", Route),
      paste("Reference:", Reference),
      paste("Media:", Media),
      paste("LOQ:", LOQ), sep = "\n")]
    bestfit.invivo[Ccompartment<0,Ccompartment:=0]

    # Determine y-axis range:
    min.conc <- 10^floor(log10(min(unique(as.numeric(this.compound.data$Value)),
                                   na.rm=T)))
    if (is.finite(mean(as.numeric(this.compound.data$LOQ)) &
                  !is.na(mean(as.numeric(this.compound.data$LOQ)))))
      min.conc <- min(min.conc,
                      mean(as.numeric(this.compound.data$LOQ))/5)
    max.conc <- 10^ceiling(log10(max(unique(as.numeric(this.compound.data$Value)),
                                     na.rm=T)))

    plotted.data <- sort(unique(this.compound.data[,(color.factor)]))
    if (plot.httk.pred)
    {
      names(colorvals)[2:(1+length(plotted.data))] <- plotted.data
    } else {
      names(colorvals) <- plotted.data
    }

if (any(!is.finite(this.graph.fits$AIC)))
{
  p <- ggplot2::ggplot(data=this.compound.data) +
    ggtitle(paste(this.compound," (",model,"): Optimizer Failed, No Curve Fit",sep=""))+
    facet_wrap(vars(Dose.nominal.units.type)) + ### added ref, media, loq
    geom_point(aes(x=Time, y=Value, color=color.factor, shape=L.O.Q.)) +
    scale_shape_manual(values = c(16, 1)) +
    # facet_wrap(facets=~vars(plot.split.factor), scales="free") +
    scale_x_continuous(limits = c(min.time,max.time)) +
    scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) +
    scale_color_manual(values = colorvals) +
    xlab("Time (h)") + ylab("Concentration (mg/L)") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    guides(fill = guide_legend(override.aes = list(color = NA)),
           color = FALSE)

} else {

  p <- ggplot2::ggplot(data=this.compound.data) +
    ggtitle(paste(this.compound," (",model,")",sep=""))+
    facet_wrap(vars(Dose.nominal.units.type)) + ### added ref, media, loq
    geom_line(data=bestfit.invivo,aes(x=time,
                                      y=Ccompartment,
                                      color=color.factor)) +
    geom_point(data=this.compound.data,
               aes(x=Time, y=Value,
                   color=color.factor,
                   shape=L.O.Q.)) +
    scale_shape_manual(values = c(16, 1)) +
    #        facet_wrap(facets=~vars(plot.split.factor), scales="free") +
    scale_x_continuous(limits = c(min.time,max.time)) +
    scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) +
    scale_color_manual(values=colorvals) +
    xlab("Time (h)") + ylab("Concentration (mg/L)") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    guides(fill = guide_legend(override.aes = list(color = NA)),
           color = FALSE)

}
    # if (any(is.finite(this.graph.fits$AIC)))
    # {
    #   p <- ggplot2::ggplot(data=this.compound.data) +
    #     ggtitle(paste(this.compound," (",model,")",sep=""))+
    #     facet_wrap(vars(Dose.nominal.units.type)) + ### added ref, media, loq
    #     geom_line(data=bestfit.invivo,aes(x=time,
    #                                       y=Ccompartment,
    #                                       color=color.factor)) +
    #     geom_point(data=this.compound.data,
    #                aes(x=Time, y=Value,
    #                    color=color.factor,
    #                    shape=shape.factor)) +
    #     scale_shape_manual(values = c(16, 1)) +
    #     #        facet_wrap(facets=~vars(plot.split.factor), scales="free") +
    #     scale_x_continuous(limits = c(min.time,max.time)) +
    #     scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) +
        # scale_color_manual(values=colorvals) +
    #     xlab("Time (h)") + ylab("Concentration (mg/L)") +
    #     theme_bw() +
    #     theme(aspect.ratio = 1)
    # } else {
    #   p <- ggplot2::ggplot(data=this.compound.data) +
    #     ggtitle(paste(this.compound," (",model,"): Optimizer Failed, No Curve Fit",sep=""))+
    #     facet_wrap(vars(Dose.nominal.units.type)) + ### added ref, media, loq
    #     geom_point(aes(x=Time, y=Value, color=color.factor, shape=shape.factor)) +
    #     scale_shape_manual(values = c(16, 1)) +
    #     # facet_wrap(facets=~vars(plot.split.factor), scales="free") +
    #     scale_x_continuous(limits = c(min.time,max.time)) +
    #     scale_y_log10(label=scientific_10,limits=c(min.conc,max.conc)) +
    #     scale_color_manual(values=colorvals) +
    #     xlab("Time (h)") + ylab("Concentration (mg/L)") +
    #     theme_bw() +
    #     theme(aspect.ratio = 1)
    # }
    if (plot.loq)
    {

      ### fixed data.table notation
      ### need to make LOQ line unique to each facet
      p <- p + geom_hline(aes(yintercept = 2 * this.compound.data[, LOQ]),
                          color = "Blue", linetype = "dashed") +
        geom_text(aes((max.time + min.time) / 2,
                      2 * this.compound.data[, LOQ],
                      label = "L.O.Q."))

      # p <- p + geom_hline(aes(yintercept=2*mean(this.compound.data[,"LOQ"])),
      #                     color="Blue",linetype="dashed") +
      #   annotate("text",x=(max.time+min.time)/2,
      #            y=2*mean(data.set[Compound==this.compound,"LOQ"]),
      #            label="Limit of Quantitiation")
    }

    print(p)
  }

  dev.off()
  return(0)
}
