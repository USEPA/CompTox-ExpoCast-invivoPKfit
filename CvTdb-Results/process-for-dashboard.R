#Clear memory
rm(list=ls())

# Selecte optimizer method based on Gil's paper:
PREFERRED.METHOD <- "L-BFGS-B"

# File with the parameter estimates by chemical and species:
PARAM.FILE <- "CvTdb_fullCoefsAIC_2023.10.30.csv"

# File with the PMIDs for each study:
CVT.FILE <-"CvTdb_fullData_20230906.csv"

# File that identifies the refernence for each PMID:
REF.FILE <- "PMID-to-REF.txt"

# Working directory:
setwd("C:/Users/jwambaug/git/invivopkfit/CvTdb-Results")

# Load the per chemical/species paramerer estimates:
param.table <- read.csv(PARAM.FILE)

# Drop bad fits:
param.table <- subset(param.table, !(model %in% "model_flat"))

# Identify the unique optimize
optim.methods <- unique(param.table$method)

# Preferentially select fits from preferred method:
param.table2 <- subset(param.table, method %in% PREFERRED.METHOD)

param.table.other <- subset(param.table, !(method %in% PREFERRED.METHOD))

for (this.chem in unique(param.table.other$Chemical))
  for (this.species in unique(param.table.other$Species))
# Do we have values for a species-chemical combination from the other method:
       if (any(param.table.other$Chemical == this.chem &
           param.table.other$Species == this.species))
# Do we NOT have values for that same combination with the preferred method:
           if (!any(param.table2$Chemical == this.chem &
               param.table2$Species == this.species))
           {
             this.subset <- subset(param.table.other, Chemical == this.chem &
                                   Species == this.species)
                                   
             # If there are multiple methods pick one at random:
             if (dim(this.subset)[1]>1) this.subset <- this.subset[
                 sample(1:dim(this.subset)[1],1), ]
             
             # Add the results from the other method:
             param.table2 <- rbind(param.table2, this.subset)
           }

# Calculate Vdist [L/kg] for two compartment
param.table2[param.table2$model %in% "model_2comp",
             "Vdist"] <- 
             param.table2[param.table2$model %in% "model_2comp", "V1"] * (1 + 
             param.table2[param.table2$model %in% "model_2comp", "k12"] / 
             param.table2[param.table2$model %in% "model_2comp", "k21"])
             
# Set reasonable sig figs:
param.table2$Vdist <- signif(param.table2$Vdist, 4)
# Calculate half-life [h]:
param.table2$Thalf <- signif(log(2)/param.table2$kelim, 4)
# Calculate clearance [L/kg/h]:
param.table2$CLtot <- signif(param.table2$Vdist * param.table2$kelim, 4)
# Calculate Css 1 mg/kg/day dose [mg/L]:
param.table2$Css <- signif(1/(24*param.table2$CLtot), 4)
# Set reasonable sig figs: 
param.table$Fgutabs <- signif(param.table$Fgutabs, 4)
 
# File emailed from PubMed when searched with all PMID's:       
pmid <- read.table(REF.FILE,sep="\t")
pmid$PMID <- unlist(lapply(strsplit(unlist(strsplit(unlist(lapply(strsplit(pmid[,1],"PMID: "),function(x) x[2])),"\\.")),";"),function(x) x[1])) 
pmid$First.Author <- lapply(strsplit(pmid[,1]," "),function(x) x[2])
for (this.month in c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
{
  which.rows <- regexpr(paste0(" ",this.month),pmid[,1])!=-1 
  pmid[which.rows,"Year.Start"] <-
    regexpr(paste0(" ",this.month),pmid[which.rows,1])-4
}
pmid$Year <- substring(pmid[,1],pmid$Year.Start,pmid$Year.Start+3)
pmid$Ref <- paste(pmid$First.Author," et al. (",pmid$Year,")",sep="")
pmid$DOI <- lapply(strsplit(unlist(lapply(strsplit(pmid[,1],"doi: "),function(x) x[2])),"\\. "),function(x) x[1])

cvt <- read.csv(CVT.FILE)

for (this.species in unique(param.table2$Species))
  for (this.chem in unique(param.table2$Chemical))
    if(any(param.table2$Species == this.species &
           param.table2$Chemical == this.chem))
    {
      NTP.data <- FALSE
      
      these.pmids <- unique(subset(cvt,subjects.species==this.species &
                                   chemicals_analyzed.dsstox_substance_id==this.chem)$documents_extraction.pmid)
      if (any(is.na(these.pmids))) these.refs <- "NTP"
      else these.refs <- NULL 
           
      if (any(!is.na(these.pmids)))
      {
        these.pmids <- these.pmids[!is.na(these.pmids)]
        these.refs <- c(these.refs,
                        unique(pmid[pmid$PMID %in% these.pmids, "Ref"]))
      } 
      param.table2[param.table2$Species == this.species &
                   param.table2$Chemical == this.chem,
                   "Ref"] <- paste(these.refs,collapse=", ")
      print(paste(these.refs,collapse=", "))
    }
param.table2[param.table2$Ref=="", "Ref"] <- "National Toxicology Program"

    
#Set reasonable sigfigs:
for (this.col in colnames(param.table2))
  if (any(!(is.na(as.numeric(param.table2[,this.col])))))
    param.table2[,this.col] <- signif(as.numeric(param.table2[,this.col]),2)
    
write.csv(param.table2, 
          file="invivoPKfit-params.for.dashboard.txt",
          row.names=FALSE)    