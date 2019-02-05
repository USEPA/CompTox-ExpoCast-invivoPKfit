## ---- include=FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = '#>')

## ----load_packages, eval = FALSE-----------------------------------------
#  library(invivoPKfit)

## ----initialize, eval = FALSE--------------------------------------------
#  TeachingDemos::char2seed("Caroline Ring")

## ----noncomp, eval=FALSE-------------------------------------------------
#  system.time(PK.fit.table.noncomp <- invivoPKfit::fit_all(data.set=pkdataset_nheerlcleaned, model="noncompartment"))
#  
#  saveRDS(PK.fit.table.noncomp , paste("output/PK_fit_table_noncomp-",Sys.Date(),".rda",sep=""))

## ----twocomp, eval=FALSE-------------------------------------------------
#  system.time(PK.fit.table.2comp  <- invivoPKfit::fit_all(pkdataset_nheerlcleaned, model="2compartment", modelfun="analytic"))
#  
#  saveRDS(PK.fit.table.2comp, paste("output/PK_fit_table_2comp-",Sys.Date(),".rda",sep=""))

## ----twcomp_plot, eval=FALSE---------------------------------------------
#  junk <- plot_conctime(PK.fit.table=PK.fit.table.2comp,
#                        data.set=pkdataset_nheerlcleaned,
#                        model="2compartment")

## ----onecomp, eval=FALSE-------------------------------------------------
#  system.time(PK.fit.table.1comp <- invivoPKfit::fit_all(data.set=pkdataset_nheerlcleaned, model="1compartment", modelfun="analytic"))
#  
#  saveRDS(PK.fit.table.1comp, paste("output/PK_fit_table_1comp-",Sys.Date(),".rda",sep=""))

## ----onecomp_plot,eval=FALSE---------------------------------------------
#  junk <- plot_conctime(PK.fit.table=PK.fit.table.1comp,
#                        data.set=pkdataset_nheerlcleaned,
#                        model="1compartment")

## ----dataoutput,eval=FALSE-----------------------------------------------
#  write.csv(pkdataset_nheerlcleaned[order(pkdataset_nheerlcleaned$Compound),],"SupTable1.txt",row.names=F)
#  write.csv(pkdataset_nheerlcleaned[order(pkdataset_nheerlcleaned$Compound),],file=paste("InVivoData-",Sys.Date(),".txt",sep=""),row.names=F)
#  save(pkdataset_nheerlcleaned,file=paste("PKdata-",Sys.Date(),".RData",sep=""))

