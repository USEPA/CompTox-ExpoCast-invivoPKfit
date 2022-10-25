#' Merge tables from fit_all from different models and select best fitting model
#'
#' This takes multiple tables of pharmcokinetic parameters (one table per model)
#' genetated by \code{\link{fit_all}}
#' and merges them into a table containing one row per chemical-species 
#' combination. Each row contains all the model
#' parameters estimated for that model, and a set of optimal parameters as 
#' identified from the AIC. This function assigns each chemical a column 
#' "Model" containing the name of he winning (lowest AIC) model and, from that
#' winning model, reports "Vdist", "kelim", "thalf", "kgutabs", "Fbioavail".
#'
#' @param fit_list a named list of PK parameter tables from \code{\link{fit_all}}
#' -- the names should be the names of the models
#'
#' @return A data.frame with one row per chemical-species combination
#'
#' @author John Wambaugh
#'
#' @export merge_model_fits 
function <- merge_model_fits(fit.list)
{
}
  fittable <- data.frame(
        DTXSID=NA,
        Compound=NA,
        CAS=NA,
        Species=NA,
        Reference=NA,
        AIC.best=NA,
        Model=NA,
        Vdist=NA,
        kelim=NA,
        kgutabs=NA,
        Fgutabs=NA)

  for (this.table.name in sort(unique(names(fit.list))))
  {
    this.table <- fit.list[[this.table.name]]
    this.table$Css <- 1/this.table$CLtot
    
    # Check for pathological fits, set AIC to Inf:
    if ("k12" %in% colnames(this.table))
    {
      # Negative k12 indicates bad fit:
      badfits <- this.table$k12
      badfits[is.na(badfits)]<-0
      badfits <- ifelse(badfits<0,TRUE,FALSE)
      this.table[badfits,"AIC"] <- Inf
    }    
   
    for (this.dtxsid in this.table$DTXSID) 
    { 
      this.compound <- subset(this.table, DTXSID==this.dtxsid)$PREFERRED_NAME[1]
      this.cas <- subset(this.table, DTXSID==this.dtxsid)$CASRN[1]

      this.subset <- subset(this.table,DTXSID == this.dtxsid)
      
      for (this.species in unique(this.subset$Species))
      {
        this.species.subset < subset(this.table, )
      if (this.dtxsid %in% fittable$dtxsid)
      {
        this.index <- which(fittable$dtxsid == this.dtxsid)
      } else this.index <- dim(fittable)[1]+1
      
      this.row <- data.frame(
        DTXSID=this.dtxsid,
        Compound=this.compound,
        CAS=this.cas,
        Species=this.species,
        Reference=NA,
        AIC.best=NA,
        Model=NA,
        Vdist=NA,
        kelim=NA,
        kgutabs=NA,
        Fgutabs=NA)

    
  for (this.species in unique(tolower(c(
    this.data.1comp$Species,
    this.data.2comp$Species,
    this.data.flat$Species))))
  {
    Vdist.1comp <- NA
    kelim.1comp <- NA
    Vdist.2comp <- NA
    kelim.2comp <- NA
    
    this.data.1comp.species <- subset(this.data.1comp,tolower(Species)==this.species)
    this.data.2comp.species <- subset(this.data.2comp,tolower(Species)==this.species)
    this.data.flat.species <- subset(this.data.flat,tolower(Species)==this.species)
    
    
    if (this.cas %in% fitsflat$CASRN)
    {
      this.data.flat.species <- subset(this.data.flat.species,
        param.value.type=="Fitted geometric mean")

    # If more than one source, require that the joint analysis worked:
     if (dim(this.data.flat.species)[1]>1)
      {
        this.data.flat.species<- subset(this.data.flat.species, 
                                         Data.Analyzed=="Joint Analysis")[1,]
      }
      if (dim(this.data.1comp.species)[1]>0) 
      {
        this.row$AIC.flat <- this.data.flat.species$AIC
        this.row$Reference <- this.data.flat.species$Reference
      }
    }
    
    if (this.cas %in% fits1comp$CASRN)
    {
      this.data.1comp.species <- subset(this.data.1comp.species,
        param.value.type=="Fitted geometric mean")

    # If more than one source, require that the joint analysis worked:
     if (dim(this.data.1comp.species)[1]>1)
      {
        this.data.1comp.species<- subset(this.data.1comp.species, 
                                         Data.Analyzed=="Joint Analysis")[1,]
      }
      if (dim(this.data.1comp.species)[1]>0) 
      {
        this.row$AIC.1comp <- this.data.1comp.species$AIC
        this.row$Vdist.1comp <- this.data.1comp.species$Vdist
        this.row$kelim.1comp <- this.data.1comp.species$kelim
        this.row$kgutabs.1comp <- this.data.1comp.species$kgutabs
        this.row$Fgutabs.1comp <- this.data.1comp.species$Fgutabs
        this.row$Reference <- this.data.1comp.species$Reference
      }
    }
    
    if (this.cas %in% fits2comp$CASRN)
    {
      this.data.2comp.species <- subset(this.data.2comp.species,
        param.value.type=="Fitted geometric mean")

    # If more than one source, require that the joint analysis worked:
      if (dim(this.data.2comp.species)[1]>1)
      {
        this.data.2comp.species <- subset(this.data.2comp.species, 
                            Data.Analyzed=="Joint Analysis")[1,]
      }
      if (dim(this.data.2comp.species)[1]>0) 
      {
        this.row$AIC.2comp <- this.data.2comp.species$AIC
        this.row$Vdist.2comp <- this.data.2comp.species$Vss
        this.row$kelim.2comp <- this.data.2comp.species$kelim
        this.row$V1.2comp <- this.data.2comp.species$V1
        this.row$k12.2comp <- this.data.2comp.species$k12
        this.row$k21.2comp <- this.data.2comp.species$k21
        this.row$kgutabs.2comp <- this.data.2comp.species$kgutabs
        this.row$Fgutabs.2comp <- this.data.2comp.species$Fgutabs
        this.row$Ralphatokelim.2comp <- this.data.2comp.species$Ralphatokelim
        this.row$Fbetaofalpha.2comp <- this.data.2comp.species$Fbetaofalpha
        this.row$alpha.2comp <- this.data.2comp.species$alpha
        this.row$beta.2comp <- this.data.2comp.species$beta
        this.row$Reference <- this.data.2comp.species$Reference
      }
    }

    
    
    if (is.na(this.row$AIC.flat)) this.row$AIC.flat <- Inf
    if (is.na(this.row$AIC.1comp)) this.row$AIC.1comp <- Inf
    if (is.na(this.row$AIC.2comp)) this.row$AIC.2comp <- Inf
    
    AICs <- c(this.row$AIC.flat,this.row$AIC.1comp,this.row$AIC.2comp)
    if (all(AICs==Inf))
    {
      this.row$Model <- "None"
      this.row$AIC.best <- NA
    } else if (this.row$AIC.2comp == min(AICs))
    {
      this.row$Model <- "2Comp"
      this.row$AIC.best <- this.row$AIC.2comp
    } else if (this.row$AIC.1comp == min(AICs))
    {
      this.row$Model <- "1Comp"
      this.row$AIC.best <- this.row$AIC.1comp
    } else {
      this.row$Model <- "Flat"
      this.row$AIC.best <- this.row$AIC.flat
    }
    
    if (this.row$Model != "Flat")
      if (this.row$Model == "1Comp")
      {
        this.row$Vdist <- this.row$Vdist.1comp
        this.row$kelim <- this.row$kelim.1comp
        this.row$kgutabs <- this.row$kgutabs.1comp
        this.row$Fgutabs <- this.row$Fgutabs.1comp
      } else 
      {
        this.row$Vdist <- this.row$Vdist.2comp
        this.row$kelim <- this.row$kelim.2comp
        this.row$kgutabs <- this.row$kgutabs.2comp
        this.row$Fgutabs <- this.row$Fgutabs.2comp
      }      
    
    fittable <- rbind(fittable,this.row)
  }
}
fittable$halflife <- log(2)/fittable$kelim

#Lets only keep chemicals where we can fit a PK model:
fittable <- subset(fittable, !is.na(Model))

```

# Set the signfifant digits to something reasonable:
```{r subset_cvt_data, eval = TRUE}
CvT.data$Value <- signif(CvT.data$Value,3)
CvT.data$Time <- signif(CvT.data$Time,3)
CvT.data$calc_loq <- signif(CvT.data$calc_loq,3)

for (this.col in c(
  "AIC.flat", "AIC.1comp", "AIC.2comp", "AIC.best", "Vdist", "kelim", "kgutabs",
  "Fgutabs", "kgutabs.1comp","kgutabs.2comp",    
  "Vdist.1comp", "kelim.1comp", "Fgutabs.1comp", "Vdist.2comp", "V1.2comp",
  "k12.2comp", "k21.2comp", "kelim.2comp", "Fgutabs.2comp", 
  "Ralphatokelim.2comp","Fbetaofalpha.2comp","alpha.2comp","beta.2comp",
  "halflife"))
  fittable[,this.col] <- signif(fittable[,this.col], 3)
```

# We have TK parameters estimated for 262 chemicals:
```{r load_physchem, eval = TRUE}
dim(fittable)[1]
length(unique(fittable$DTXSID))

write.csv(fittable,file="SupTable-TKFits.txt",row.names=FALSE)
```
