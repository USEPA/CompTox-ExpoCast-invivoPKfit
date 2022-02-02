load("TKparams-2020-02-28.RData")

max(all_study_1c_fits.mean$kelim)
min(all_study_1c_fits.mean$kelim)
min(all_study_1c_fits.mean$Vdist)
max(all_study_1c_fits.mean$Vdist)

all_study_1c_fits.better <- subset(all_study_1c_fits.mean,
  kelim < 1000 & Vdist > .1 & Vdist < 1000)
  
all_compound_1c_fits <- NULL
for (this.cas in unique(all_study_1c_fits.better$CAS))
{
  this.subset <- subset(all_study_1c_fits.better, CAS==this.cas)
  if (dim(this.subset)[1]==1) all_compound_1c_fits <- 
    rbind(all_compound_1c_fits, this.subset)
  else {
    if (any(regexpr(",",this.subset$Reference)!=-1))
    {
      all_compound_1c_fits <- 
        rbind(all_compound_1c_fits, this.subset[
          regexpr(",",this.subset$Reference)!=-1,])
    } else {
      new.row <- this.subset[1,]
      new.row$Reference <- paste("Mean of",paste(
        this.subset$Reference,collapse=", "))
      new.row$Data.Analyzed <- paste(
        this.subset$Reference,collapse=", ")
      for (this.col in 
        c("kelim",
          "Vdist",
          "Fgutabs",
          "kgutabs",
          "CLtot",
          "Css",
          "halflife",
          "tpeak.oral",
          "Cpeak.oral.1mgkg")) new.row[,this.col] <- 
        mean(this.subset[,this.col],na.rm=T)
      for (this.col in 
        c("LogLikelihood",
          "AIC")) new.row[,this.col] <- NA 
      all_compound_1c_fits <- rbind(all_compound_1c_fits, new.row)       
    }
  }
}

save(all_compound_1c_fits,all_study_1c_fits.better,
  file="CvT-invivoPKfit-102220.RData")