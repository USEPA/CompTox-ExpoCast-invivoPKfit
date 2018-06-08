library(httk)
library(gdata)
library(data.table)
library(stringr)

#Set perl location for later use by gdata::read.xls()
#perldir <- 'C:/Users/CRing/Applications/Perl/bin/perl.exe'

#Pull in the TNO data, included in httk package
#Possible dose types
#mg/kg po: oral administration
#mg/kg iv: intravenous administration
#mg/kg sc: subcutaneous administration
#5x mg/kg sc: 5 times recommended dose, subcutaneous administration
#Note: For some reason all of the rows of chem.invivo.PK.data are duplicated??
#Remove duplicates.
#TNO.dat <- unique(as.data.table(httk::chem.invivo.PK.data))
#diclofenac data is duplicated for some reason. Remove duplicates.
#TNO.dat <- subset(TNO.dat, CAS!="64118-84-9")
#load("L:/Lab/NCCT_ExpoCast/ExpoCast2016/HTTKDataTable")
chem.invivo.PK.data[chem.invivo.PK.data[,"Media"]=="Plasma conc","Media"]<-"Plasma concentration"
TNO.dat <- as.data.table(chem.invivo.PK.data)

#Note that all of Species should be "rat" but at least one entry has a trailing space.
#remove the trailing space.
TNO.dat[, Species:=gsub(x=Species,pattern=" ",replacement="")]

#Strip spaces from column names
setnames(TNO.dat, names(TNO.dat), gsub(x=names(TNO.dat),
                                         pattern=' ',
                                         replacement='.'))

#Exclude anything sc because neither po nor iv
#This includes caffeine and diltiazem
TNO.dat <- subset(TNO.dat, !grepl(x=Dose.Units.and.Type,
                                 pattern="sc"))

#Set dose units
TNO.dat[, Units:='mg/kg']
#Get dose routes
TNO.dat[, Route:=gsub(x=Dose.Units.and.Type,
                       pattern='mg/kg ',
                       replacement='',
                       fixed=TRUE)]

#Take only the chemicals for which the 1-compartment model can be parameterized
#in humans or in rats
TNO.dat <- subset(TNO.dat,CAS %in% unique(c(httk::get_cheminfo(model="1compartment",
                                                               species="Human"),
                                              httk::get_cheminfo(model="1compartment",
                                                           species="Rat"))))
#Take only the rows for which "Value" (the measured concentration) was greater
#than zero
TNO.dat <- subset(TNO.dat,Value>0)

#Add a column to the data set denoting the source of these data:
#here, the TNO Lit. Review
TNO.dat[, Source:='TNO Lit. Review']

# #Change time points of 0 to 1/10 of the time points following them, for easier
# #fitting
# TNO.dat[which(TNO.dat$Time==0),
#          Time:=TNO.dat[which(TNO.dat$Time==0)+1,Time]/10]
# Not sure if TNO pulled MEHP or DEHP data
#so exclude anything labeled Mono ethylhexyl phthalate
TNO.dat <- subset(TNO.dat,Compound!="Mono ethylhexyl phthalate")


# Two weights are reported for Tokuma 1988, so average them (we get an NA otherwise):
TNO.dat[Reference=="Tokuma 1988",Species.Weight:=(0.24+0.18)/2]

#Strip parentheses from TNO time units and media units
TNO.dat[, Time.Units:="h"]
TNO.dat[, Media.Units:="ug/mL"]

# Add a fake Limit of Quantitation (LOQ) for literature data such that all 
# points are 2x above it:
TNO.dat[,LOQ:=0.45*min(Value),by="CAS,Reference"]

#Now read in the RTI data

#Read in the first sheet, "Test Articles"
#This contains data on each chemical: its name, CAS, MW (g/mol),
#the supplier from which it was obtained,
#its lot number, and manufacturer's purity
RTI.chems <- gdata::read.xls("RTI/RTIdata.xlsx",
                      stringsAsFactors=F,
                      sheet='Test Articles')
#Exclude any chemicals with a missing CAS
RTI.chems <- subset(RTI.chems,CAS..!="")
#Convert to data.table
RTI.chems <- as.data.table(RTI.chems)

#Read in sheet 4 of the RTI data, "Dose Administration"
#Contains details on dosing of each subject
#Subject ID, dose route, which compound was dosed,
#the weight of the dose formulation delivered (g),
#grams of the actual compound that were dosed,
#and the actual dose rate in mg/kg
RTI.dose <- as.data.table(gdata::read.xls("RTI/RTIdata.xlsx",
                     stringsAsFactors=F,
                     sheet='Dose Administration'))
#Strip slashes from column names for usability later (replace with .)
setnames(RTI.dose, names(RTI.dose),
         gsub(x=names(RTI.dose),
              pattern='/',
              replacement='.',
              fixed=TRUE))
#Fix two compound names
RTI.dose[Test.Article=="Pyrithiobac Na",
         Test.Article:="Pyrithiobac sodium"]
RTI.dose[Test.Article=="Diazoxon",
         Test.Article:="Diazinon-o-analog"]

#Now for each chemical in RTI.data,
#read in the data from the appropriate sheet of the Excel spreadsheet.

RTI.readfun <- function(compound){
  if (compound=='Diazinon-o-analog') {
    #The sheet is named something different from the compound
    compound.tmp <- 'Diazoxon Metabolite IMP'
  } else if (compound=='Pyrithiobac sodium'){
    compound.tmp <- 'Pyrithiobac Na'
  } else{
    compound.tmp <- compound
  }

  tryCatch(
  this.data <- as.data.table(gdata::read.xls("RTI/RTIdata.xlsx",
                        stringsAsFactors=F,
                        header=FALSE,
                        sheet=compound.tmp,
                        skip=1)),
  error = function(err){
    if (err$message == "invalid 'file' argument"){
      #try adding a trailing space to the sheet name
      tryCatch(this.data <- as.data.table(gdata::read.xls("RTI/RTIdata.xlsx",
                                          stringsAsFactors=F,
                                          header=FALSE,
                                          sheet=paste0(compound.tmp, ' '),
                                          skip=1)), #skip the first line because it's just a title,
               error = function(err){
                 stop(err) #if that doesn't work, rethrow error
               }
      )
    }
  }
  )

  #The data will be in a 8-column table.
  #Column 1 is time.
  #Row 1 specifies route.
  #Row 2 gives subject IDs.
  #Rows 3-end are actual data.
  #Need to reshape this.
  setnames(this.data, 'V1', 'Time')
  # Find the limit of quantitation (LOQ):
  LOQ <- this.data[regexpr("BLOQ",this.data$Time)!=-1,"Time",with=F]
  # Extract the number:
  LOQ <- as.numeric(strsplit(strsplit(as.character(LOQ),"\\(")[[1]][2]," ")[[1]][1])
  # Convert to mg/L:
  LOQ <- LOQ/1000
  
  this.data[, Time:=as.numeric(Time)]
  #the 'Time (h)' label in row 3 will become NA
  #Also any notes at the bottom, which will have gone into V1,
  #will become NA.

  #Now take each column of time/value data in rows 3:end.
  #Add a column for Route and Subject from rows 1 and 2.
  #then rbind them together.
  this.data.shape <- rbindlist(lapply(paste0('V',2:7),
         function(x) data.table(Time=this.data[3:nrow(this.data), Time],
                           Value=as.numeric(this.data[3:nrow(this.data),
                                     eval(parse(text=x))])/1000, #convert ng/mL to mg/L
                           Route=this.data[1, eval(parse(text=x))],
                           Subject=this.data[2, eval(parse(text=x))])))

  #Remove any rows where Time is NA; it means there's no data in that row
  this.data.shape <- subset(this.data.shape, !is.na(Time))

  #Add a column with the compound name
  this.data.shape[, Test.Article:=compound]
  
  #Add a column with the LOQ:
  this.data.shape[,LOQ:=LOQ]
  
  return(this.data.shape)
}

#Now read in data for each chemical and rbind into one big data table
RTI.data <- rbindlist(lapply(unique(RTI.dose$Test.Article),
                   RTI.readfun))

#Diazoxon metabolite IMP is formated different for LOQ:
RTI.data[is.na(RTI.data$LOQ),LOQ:=0.1]

#Get species weight and dose by matching with RTI.dose
#for concordance of subject IDs, change 00 to 0 in RTI.dose
RTI.dose[, Subject:= gsub(x=Subject,
                          pattern='00',
                          replacement='0')]
RTI.dose[Route=="Gavage", Route:='PO']
#Now do the merge
RTI.data <- merge(RTI.data,
                  RTI.dose[, .(Subject,
                               Route,
                               Test.Article,
                               Subject.Wt..g.,
                               Actual.Dose.Rate..mg.kg.)],
                  by=c('Subject', 'Route', 'Test.Article'))

RTI.data <- merge(RTI.data,
                  RTI.chems[, .(Test.Article, CAS..)],
                  by='Test.Article')

setnames(RTI.data, c('Test.Article',
                     'Subject.Wt..g.',
                     'Actual.Dose.Rate..mg.kg.',
                     'CAS..'),
         c('Compound',
           'Species.Weight',
           'Dose',
           'CAS'))
RTI.data[,Species.Weight:=Species.Weight/1000]
RTI.data[,Reference:="RTI 2015"]
RTI.data[,Source:="RTI 2015"]
RTI.data[,Species:="rat"]
RTI.data[,Media:="Plasma concentration"]
RTI.data[,Time.Units:="h"]
RTI.data[,Species.Weight.Units:="kg"]
RTI.data[,Media.Units:="ug/mL"]
RTI.data[, Route:=tolower(Route)]
RTI.data[, Dose.Units.and.Type:=paste('mg/kg', Route)]
RTI.data[,Units:="mg/kg"]

#Now the data set contains both TNO and RTI data.
#Next: Add NHEERL data.

#First, read in NHEERL list of chemicals and doses.
NHEERL.chems <- as.data.table(read.xls("NHEERL/ToxCast Bioavailability List of chemicals and doses_021915.xlsx",
                         stringsAsFactors=F,
                         skip=2))

NHEERL.LOQs <- read.xls("NHEERL/2014-008 LOQs.xlsx",
                         stringsAsFactors=F,
                         skip=2)

#Read in table of compound names and CAS numbers for NHEERL studies
CAS.table <- as.data.table(read.xls("evaluation-chems-summary-072214.xlsx",
                      stringsAsFactors=F))



nheerl_data_parse <- function(nheerl.data,
                              NHEERL.chems,
                              CAS.table)
{
  nheerl.data <- copy(nheerl.data) #so it acts as though passed by value, rather than by reference
  setnames(nheerl.data,
           names(nheerl.data),
           c('info', 'Value'))
  
  #Harmonize chemical names with those in the NHEERL.chems and CAS.table
  #replace HCL, HCl with hydrochloride
  nheerl.data$info <- gsub(x=nheerl.data$info,
                           pattern='HCL|HCl',
                           replacement='hydrochloride')
  
  #Replace PFOA with its full name (AKA perfluorooctanoic acid)
  nheerl.data$info <- gsub(x=nheerl.data$info,
                           pattern='PFOA',
                           replacement='Pentadecafluorooctanoic acid')
  #And harmonize the name in the CAS table
  CAS.table[Compound=='Perfluorooctanoic acid',
            Compound:='Pentadecafluorooctanoic acid']
  
  #Bisphenol-A = Bisphenol A
  CAS.table[Compound=='Bisphenol-A',
            Compound:='Bisphenol A']
  
  #Bensulide hydrochloride == bensulide
  nheerl.data$info <- gsub(x=nheerl.data$info,
                           pattern='Bensulide hydrochloride',
                           replacement='Bensulide')
  
  nheerl.data$info <- gsub(x=nheerl.data$info,
                           pattern='S-bioallethrin',
                           replacement='S-Bioallethrin')
  
  #Somebody misspelled formetanate in the NHEERL.chems file
  #so fix it
  NHEERL.chems$Chemical <- gsub(x=NHEERL.chems$Chemical,
                           pattern='Formetane',
                           replacement='Formetanate')
  #They also did it for one subject in the nheerl.data file
  #so fix it there too
  nheerl.data$info <- gsub(x=nheerl.data$info,
                                pattern='Formetane',
                                replacement='Formetanate')
  
  #Replace two spaces with one space
  nheerl.data$info <- gsub(nheerl.data$info,
                           pattern='  ',
                           replacement=' ')
  
  #Pattern matching in nheerl.data$info
  #Pattern goes like this:
  #[Chemical name] [Subject ID] [Route] [Time] min
  #Chemical name may have spaces in it, so can't just split on spaces.
  #Need to match for any of the list of NHEERL.chems, then split the rest on spaces.
  #Paste NHEERL chems into a pipe delimited bunch.
  nheerlchems.pipe <- paste(NHEERL.chems$Chemical, collapse='|')
  #Then do regexpr and get regmatches
  nheerl.data[, Compound:=regmatches(nheerl.data$info,
                                     regexpr(text=nheerl.data$info,
                               pattern=nheerlchems.pipe,
                               perl=TRUE))]
  #Remove the compound names from nheerl.data$info
  nheerl.data$info <- gsub(x=nheerl.data$info,
                           pattern=nheerlchems.pipe,
                           replacement='')
  #Now the remaining info ought to be space delimited as
  #[Subject ID] [Route] [Time]
  #Sometimes it isn't though -- spaces are missing in a few cases -- so parse for that
  #It should be pretty easy -- it'll be
  #[numbers] [optional space] [letters] [optional space] [numbers] [optional space] 'min'
  r <- '(\\d+)\\s?([a-zA-Z]+)\\s?(\\d+)\\s?min' #translated into perl
  
  tmp<-stringr::str_match_all(nheerl.data$info, r)
  #returns a list of matrices, each matrix with columns:
  #entire string with match; capture group 1; capture group 2; capture group 3
  nheerl.data[, Subject:=sapply(tmp, function(x) x[,2])]
  nheerl.data[, Time:=sapply(tmp,
                             function(x) as.numeric(x[,4]))/60] #convert time from minutes to hours
  nheerl.data[, Route:=sapply(tmp, function(x) x[,3])]
  #Assign dose units and type based on route
  nheerl.data[, Route:=tolower(Route)]
  nheerl.data[Route=='oral', Route:='po']
  nheerl.data[, Dose.Units.and.Type:=paste('mg/kg', Route)]
  
  #Get dose from route and compound name
  setnames(NHEERL.chems, 'Chemical', 'Compound')
  nheerl.data <- merge(nheerl.data,
                       NHEERL.chems[, .(Compound,
                                        Oral.Dose.mg.kg,
                                        IV.Dose.mg.kg)],
                       by='Compound')
  #Get CAS from compound
  nheerl.data <- merge(nheerl.data,
                       CAS.table[, .(Compound, CAS)],
                       by='Compound')
  
  nheerl.data[Route=='iv', Dose:=as.numeric(IV.Dose.mg.kg)]
  nheerl.data[Route=='po', Dose:=as.numeric(Oral.Dose.mg.kg)]
  nheerl.data[, IV.Dose.mg.kg:=NULL]
  nheerl.data[, Oral.Dose.mg.kg:=NULL]
  
  nheerl.data[, Value:=as.numeric(Value)/1000] #convert ng/mL to mg/L
  nheerl.data[,Reference:="NHEERL 2015"]
  nheerl.data[,Source:="NHEERL 2015"]
  nheerl.data[,Species:="rat"]
  nheerl.data[,Media:="Plasma concentration"]
  nheerl.data[,Time.Units:="h"]
  nheerl.data[,Species.Weight.Units:="kg"]
  nheerl.data[,Media.Units:="ug/mL"] #equivalent to mg/L
  nheerl.data[,Units:="mg/kg"]
  # Assume that LOQ is 99% of lowest reported value:
  nheerl.data[,LOQ:=0.99*min(Value,na.rm=T),by=Compound]
  for (chemical in NHEERL.LOQs$Analyte) nheerl.data[Compound==chemical,LOQ:=NHEERL.LOQs[NHEERL.LOQs$Analyte==chemical,"LOQ..ng.mL"]/1000]
  return(nheerl.data)
}

#Read in data manually cleaned of subjects labeled "probably misdosed"
nheerl.data.cleaned <- as.data.table(read.xls("NHEERL/nheerl_data_cleaned.xlsx",
                                              stringsAsFactors=F))
#Parse it
nheerl.data.cleaned <- nheerl_data_parse(nheerl.data.cleaned,
                                         NHEERL.chems,
                                         CAS.table)

#Do the same with the original data
#(including all subjects labeled "probably misdosed", etc.)
#I reformatted it into a single spreadsheet for ease of use
nheerl.data.orig <- as.data.table(read.xls("NHEERL/nheerl_orig_reformatted.xlsx",
                                           stringsAsFactors=F))
#Parse it
nheerl.data.orig <- nheerl_data_parse(nheerl.data.orig,
                                         NHEERL.chems,
                                         CAS.table)
data.set <- rbind(TNO.dat,
                  RTI.data,
                  nheerl.data.cleaned,
                  fill=TRUE)

pkdataset_nheerlcleaned <- data.set
#The data sets are already constructed and are part of the invivoPKfit package
#Simazine LOQ is wrong:
pkdataset_nheerlcleaned[pkdataset_nheerlcleaned$Compound == "Simazine","LOQ"]<-10^-3
#Boscalid LOQ is wrong:
pkdataset_nheerlcleaned[pkdataset_nheerlcleaned$Compound == "Boscalid","LOQ"]<-5*10^-3
#Remove triclosan oral data:
pkdataset_nheerlcleaned <- subset(pkdataset_nheerlcleaned,Route=="iv" | Compound != "Triclosan")
#Remove Bisphenol A oral data:
pkdataset_nheerlcleaned <- subset(pkdataset_nheerlcleaned,Route=="iv" | Compound != "Bisphenol A")
#Remove imidacloprid oral data:
pkdataset_nheerlcleaned <- subset(pkdataset_nheerlcleaned,Compound != "Imidacloprid")
# Add new imidacloprid oral data:
new.imidaclopridoral <- gdata::read.xls("NHEERL/RevisedImidacloprid.xlsx",
                      stringsAsFactors=F)
new.imidaclopridoral$CAS<-"138261-41-3"
new.imidaclopridoral$Reference<-"NHEERL 2015"
new.imidaclopridoral$Species<-"Rat"
new.imidaclopridoral$Species.Weight<-NA
new.imidaclopridoral$Species.Weight.Units<-"kg"
new.imidaclopridoral$Dose<-5
new.imidaclopridoral$Dose.Units.and.Type<-"mg/kg po"
new.imidaclopridoral$Time<-new.imidaclopridoral$Time/60
new.imidaclopridoral$Time.Units<-"h"
new.imidaclopridoral$Media<-"Plasma concentration"
new.imidaclopridoral$Media.Units<-"ug/mL"
new.imidaclopridoral$Units<-"mg/kg"
new.imidaclopridoral$Source<-"NHEERL 2015"
new.imidaclopridoral$Value<-as.numeric(new.imidaclopridoral$Conc)/1000
new.imidaclopridoral$LOQ<-.001
new.imidaclopridoral$info<-"Rerun 12/17"
pkdataset_nheerlcleaned <- rbind(pkdataset_nheerlcleaned,new.imidaclopridoral[,colnames(pkdataset_nheerlcleaned)])






# Make sure there is no duplicated data:
pkdataset_nheerlcleaned <- subset(pkdataset_nheerlcleaned,!duplicated(pkdataset_nheerlcleaned))
#Need subject ids to do NCA of Tokuma 1988 data:
pkdataset_nheerlcleaned[pkdataset_nheerlcleaned$Reference=="Tokuma 1988"&pkdataset_nheerlcleaned$Route=="iv","Subject"] <-  c(rep(1,7),rep(2,6))
#Need subject ids to do NCA of Tokuma 1988 data:
pkdataset_nheerlcleaned[pkdataset_nheerlcleaned$Reference=="Tokuma 1988"&pkdataset_nheerlcleaned$Route=="po","Subject"] <-  c(rep(3,7),rep(4,8))
#Doses were swapped for Bosentan:
pkdataset_nheerlcleaned[Reference=="Treiber 2004"&Dose.Units.and.Type=="mg/kg iv","Dose"] <- 1
pkdataset_nheerlcleaned[Reference=="Treiber 2004"&Dose.Units.and.Type=="mg/kg po","Dose"] <- 10
#Dose was wrong for Metoprolol:
pkdataset_nheerlcleaned[Compound=="Metoprolol","Dose"] <- 15
# Models don't handle time 0 data well:
pkdataset_nheerlcleaned<-subset(pkdataset_nheerlcleaned,Time>0)
# Remove last Binkerd 1988 time point 
pkdataset_nheerlcleaned<-subset(pkdataset_nheerlcleaned,Compound!="Valproic acid"|Time!=24)

devtools::use_data(pkdataset_nheerlcleaned, overwrite=TRUE)
#
#
#Make an alternate version with the non-cleaned NHEERL daata
#pkdataset_nheerlorig <- rbind(TNO.dat,
#                                          RTI.data,
#                                          nheerl.data.orig,
#                                          fill=TRUE)
#
#pkdataset_nheerlorig <- subset(pkdataset_nheerlorig,
#                   !is.na(Value)) #remove all NA data      

#devtools::use_data(pkdataset_nheerlorig, overwrite=TRUE)
