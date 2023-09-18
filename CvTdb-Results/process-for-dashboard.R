#Clear memory
rm(list=ls())

# Selecte optimizer method based on Gil's paper:
PREFERRED.METHOD <- "L-BFGS-B"

# Working directory:
setwd("C:/Users/jwambaug/git/invivopkfit/CvTdb-Results")

# Load the per chemical/species paramerer estimates:
param.table <- read.csv("CvTdb_fullCoefs_20230906.csv")

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
             
# Set reasonabe sig figs:
param.table2$Vdist <- signif(param.table2$Vdist, 4)
# Calculate half-life [h]:
param.table2$Thalf <- signif(log(2)/param.table2$kelim, 4)
# Calculate clearance [L/kg/h]:
param.table2$CLtot <- signif(param.table2$Vdist * param.table2$kelim, 4)
# Calculate Css 1 mg/kg/day dose [mg/L]:
param.table2$Css <- signif(1/(24*param.table2$CLtot), 4)

             