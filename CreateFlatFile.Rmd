---
title: "Create Flat File"
author: "John Wambaugh"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{TKQSPR Bakeoff}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

#Load the curve fits from invivoPKfit (Chris Cook):
```{r load_compartmental_model_parameters, eval = TRUE}
setwd("C:/Users/jwambaug/git/tkqsars/")
fitsflat <- read.csv("PK.fit.table.flat.03012022.csv")
fits1comp <- read.csv("PK.fit.table.1comp.03012022.csv")
#fits1comp <- subset(fits1comp, is.finite(AIC) & AIC < 10^4)
fits2comp <- read.csv("PK.fit.table.2comp.03012022.csv")
#fits2comp <- subset(fits2comp, is.finite(AIC) & AIC < 10^4)
# Units on Volumes:
fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"Vdist"] <-
  fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"Vdist"] /
  1000
fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"CLtot"] <-
  fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"CLtot"] /
  1000
fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"Css"] <-
  fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"Css"] *
  1000
fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"Cpeak.oral.1mgkg"] <-
  fits1comp[regexpr("geometric mean",fits1comp$param.value.type)!=-1,"Cpeak.oral.1mgkg"] *
  1000
fits2comp[regexpr("geometric mean",fits2comp$param.value.type)!=-1,"V1"] <-
  fits2comp[regexpr("geometric mean",fits2comp$param.value.type)!=-1,"V1"] /
  1000
fits2comp[regexpr("geometric mean",fits2comp$param.value.type)!=-1,"Vss"] <-
  fits2comp[regexpr("geometric mean",fits2comp$param.value.type)!=-1,"Vss"] /
  1000
fits2comp[regexpr("geometric mean",fits2comp$param.value.type)!=-1,"CLtot"] <-
  fits2comp[regexpr("geometric mean",fits2comp$param.value.type)!=-1,"CLtot"] /
  1000
fits2comp$Css <- 1/fits2comp$CLtot

# Negative k12 indicates bad fit:
badfits <- fits2comp$k12
badfits[is.na(badfits)]<-0
badfits <- ifelse(badfits<0,TRUE,FALSE)
fits2comp[badfits,"AIC"] <- Inf


#How many chemicals have 1 or 2 compartment fits:
length(unique(c(fitsflat$CAS,fits1comp$CAS,fits2comp$CAS)))
```

# Concatanate the curve fits into a single table (should be part of invivoPKfit eventually):
```{r merge_curve_fits, eval = TRUE}
fittable <- NULL
for (this.dtxsid in unique(c(fits1comp$DTXSID,fits2comp$DTXSID,fitsflat$DTXSID)))
{
  if (this.dtxsid %in% fits1comp$DTXSID) 
  { 
    this.compound <- subset(fits1comp,DTXSID==this.dtxsid)$PREFERRED_NAME[1]
    this.cas <- subset(fits1comp,DTXSID==this.dtxsid)$CASRN[1]
  } else if (this.cas %in% fits2comp$CAS) {
    this.compound <- subset(fits2comp,DTXSID==this.dtxsid)$PREFERRED_NAME[1]
    this.cas <- subset(fits2comp,DTXSID==this.dtxsid)$CASRN[1]
  } else {
    this.compound <- subset(fitsflat,DTXSID==this.dtxsid)$PREFERRED_NAME[1]
    this.cas <- subset(fitsflat,DTXSID==this.dtxsid)$CASRN[1]
  }
    
  this.data.1comp <- subset(fits1comp,DTXSID==this.dtxsid)
  this.data.2comp <- subset(fits2comp,DTXSID==this.dtxsid)
  this.data.flat <- subset(fitsflat,DTXSID==this.dtxsid)
  
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
    
    this.row <- data.frame(
      DTXSID=this.dtxsid,
      Compound=this.compound,
      CAS=this.cas,
      Species=this.species,
      Reference=NA,
      AIC.flat=NA,
      AIC.1comp=NA,
      AIC.2comp=NA,
      AIC.best=NA,
      Model=NA,
      Vdist=NA,
      kelim=NA,
      kgutabs=NA,
      Fgutabs=NA,
      Vdist.1comp=NA,
      kelim.1comp=NA,
      kgutabs.1comp=NA,
      Fgutabs.1comp=NA,
      Vdist.2comp=NA,
      V1.2comp=NA,
      k12.2comp=NA,
      k21.2comp=NA,
      kelim.2comp=NA,
      kgutabs.2comp=NA,
      Fgutabs.2comp=NA,
      Ralphatokelim.2comp=NA,
      Fbetaofalpha.2comp=NA,
      alpha.2comp=NA,
      beta.2comp=NA)
    
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
