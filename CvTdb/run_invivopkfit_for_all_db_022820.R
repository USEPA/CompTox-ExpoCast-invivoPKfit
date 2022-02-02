library(data.table)
library(invivoPKfit)
library(scales)
library(ggplot2)
library(grid)

## USEFUL PLOT FUNCTIONS:

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, heights=NULL, widths=unit(rep_len(1, cols), "null")) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    if (!is.null(heights))
    {
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout),heights=heights,widths=widths)))
    } else {
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout),widths=widths)))
    }
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

scientific_10 <- function(x) {                                  
  out <- gsub("1e", "10^", scientific_format()(x))              
  out <- gsub("\\+","",out)                                     
  out <- gsub("10\\^01","10",out)                               
  out <- parse(text=gsub("10\\^00","1",out))                    
}  

## ANALYSIS:
load("L:/Lab/NCCT_ExpoCast/SayreData/res_for_invivopkfit.Rdata")
series_res_set <- as.data.table(series_res_set)
# Set the LOQ based on chemical, paper, and species:
series_res_set[,LOQ:=min(as.numeric(Value),na.rm=T)*0.9,by=.(CAS,Reference,Species,Route)]
# We assume each source had a single analytical chemistry lab for which we want
# to estimate the standard deviation (Cvt uses reference to identify specific 
# animals):
series_res_set[,Reference:=Source]

all_study_1c_fits  <- fit_all(series_res_set,model="1compartment",modelfun = 'analytic')


all_study_1c_fits.good <- all_study_1c_fits[is.finite(AIC) & LogLikelihood!=-99999]
all_study_1c_fits.mean <- all_study_1c_fits.good[param.value.type == "Fitted geometric mean"]
fits.joint <- all_study_1c_fits.mean[Data.Analyzed == "Joint Analysis"]
fits.single <-  all_study_1c_fits.mean[!(CAS %in% fits.joint$CAS)]
fits.single <- fits.single[,lapply(.SD, function(x) ifelse(is.character(x),"Many",mean(x))),by=.(CAS,Compound,Species)]
all.fits <- rbind(fits.joint,fits.single)



## PLOTS:
 Fig1a <- ggplot(all.fits, aes(x=Vdist)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#FF6666")+ 
 xlab(expression(paste(V[dist]," (L)"))) +
  ylab("Probability Density")+
   scale_x_log10(label=scientific_10,limits=c(10^-2,10^6))+
   labs(title="Volume of Distribution")+
    theme_bw()+
    theme( text  = element_text(size=16))+
    annotate("text", x=3*10^-2,y=0.4,size=8,label="a")
  
 Fig1b <- ggplot(all.fits, aes(x=kelim)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#66FF66")+ 
 xlab(expression(paste(k[elim]," (1/h)"))) +
  ylab("Probability Density")+
   scale_x_log10(label=scientific_10,limits=c(10^-8,10^5))+
   labs(title="Elimination Rate")+
    theme_bw()+
    theme( text  = element_text(size=16))+
   annotate("text", x=10^-7,y=0.35,size=8,label="b")
   
   
 Fig1c <- ggplot(all.fits, aes(x=Fgutabs)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#6666FF")+ 
 xlab(expression(paste(F[bio]))) +
  ylab("Probability Density")+
   scale_x_log10(label=scientific_10,limits=c(10^-2,1.1))+
   labs(title="Oral Bioavailability")+
    theme_bw()+
    theme( text  = element_text(size=16))+
    annotate("text", x=1.5*10^-2,y=6,size=8,label="c")
  

    
 Fig1d <- ggplot(all.fits, aes(x=kgutabs)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#666666")+ 
 xlab(expression(paste(k[gutabs]," (1/h)"))) +
  ylab("Probability Density")+
   scale_x_log10(label=scientific_10,limits=c(10^-8,10^4))+
   labs(title="Oral Absorption Rate")+
    theme_bw()+
    theme( text  = element_text(size=16))+
    annotate("text", x=10^-7,y=0.36,size=8,label="d")
  

multiplot(Fig1a,Fig1c,Fig1b,Fig1d,cols=2,widths=c(1.75,1.75))

save.image(file=paste("TKparams-",Sys.Date(),".RData",sep=""))
write.csv(all_study_1c_fits.good,file="All-fits-joint-and-separate.txt",row.names=F)
write.csv(all.fits,file="All-fits-one-per-chemical.txt",row.names=F)

library(gridExtra)
grid.arrange(Fig1a, Fig1b, Fig1c, Fig1d, ncol=2, top = "Distributions of Estimated TK Parameters")




