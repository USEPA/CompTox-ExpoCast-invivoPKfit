get_starts_2comp <- function(data,
                             par_DF){

  kelim <- NA_real_
  kgutabs <- NA_real_
  V1 <- NA_real_
  Fgutabs_V1 <- NA_real_
  k12 <- NA_real_
  k21 <- NA_real_
  Fgutabs <- NA_real_
  Rblood2plasma <- 1

  # Get starting Concs from data

  #Work only with detects for these rough estimates
  tmpdat <- subset(data,
                   Detect %in% TRUE)

  #Split into IV and PO
  ivdat <- subset(tmpdat,
                  Route %in% "iv")
  podat <- subset(tmpdat,
                  Route %in% "oral")

  # Quick and dirty:
  #IV data estimates, if IV data exist
if(nrow(ivdat)>0){
#get "elbow"
 elbow_time <- get_elbow(x = ivdat$Time,
                         y = log10(ivdat$Conc/ivdat$Dose))
#split into late and early phases at elbow
 late_iv <- subset(ivdat, Time >= elbow_time)
 early_iv <- subset(ivdat, Time <= elbow_time)
#get slope of late phase

 #

  #assume that midpoint of time is one half-life, so kelim = 0.693/(midpoint of time).
  halflife <- mean(range(ivdat$Time))
  kelim <- 0.693/halflife

  #V1: extrapolate back from conc at min time at a slope of -kelim to get the intercept
  #then V1 = 1/intercept
  C_tmin <- with(subset(ivdat, Time == min(Time)),
                 median(Conc/Dose))
  A <- C_tmin + kelim*min(ivdat$Time)
  V1 <- 1/A
}

  if(nrow(podat)>0){
    #if PO data exist, then we can get ka and Fgutabs/V1
    #get peak time
    tCmax <- get_peak(x = podat$Time,
                      y = podat$Conc/podat$Dose)
    tmax <- tCmax[[1]]
    Cmax <- tCmax[[2]]

    #assume peak time occurs at 1 absorption halflife
    #so kgutabs = 0.693/tmax
    kgutabs <- 0.693/tmax

    #if no IV data, then calculate kelim from oral data
    if(nrow(ivdat)==0){
      #and assume that midpoint of time is one half-life, so kelim = 0.693/(midpoint of time).
      halflife <- mean(range(podat$Time))
      kelim <- 0.693/halflife
    }

    #then extrapolate back from Cmax to time 0 with slope -kelim
    Fgutabs_V1 <- (Cmax + kelim*tmax)*(kgutabs - kelim)/(kgutabs)

    if(nrow(ivdat)>0){
      #if we had IV data, then we had a V1 estimate, so we can estimate Fgutabs too
      Fgutabs <- Fgutabs_V1 * V1
    }
  }

  #arbitrary until I implement the rest of this
  k21 <- 0.1
  k12 <- 0.5

  starts <- c("kelim" = kelim,
              "kgutabs" = kgutabs,
              "k12" = k12,
              "k21" = k21,
              "V1" = V1,
              "Fgutabs_V1" = Fgutabs_V1,
              "Fgutabs" = Fgutabs,
              "Rblood2plasma" = Rblood2plasma)

par_DF$start <- starts[par_DF$param_name]

  return(par_DF)
}
