get_starts <- function(fitdata,
                       this.dtxsid,
                       modelfun,
                       model,
                       this.reference = NULL,
                       suppress.messages = FALSE,
                       sig.figs = 5){
  #Initialize starting guesses: to be refined
  if(model %in% "flat"){
    opt.params <- c("A" = mean(fitdata$Value, na.rm = TRUE))
  }else{
    #for non-flat models....
    if (model %in% "1compartment") {
      opt.params <- c("kelim" = 0.25,
                      "Vdist" = 5.56)
    }else if(model %in% "2compartment"){
      opt.params <- c("kelim" = 0.25,
                      "V1" = 1,
                      "Ralphatokelim" = 1.2,
                      "Fbetaofalpha" = 0.25)
    }

    if('po' %in% fitdata$Route){
      #if we have oral data
      #then add Fgutabs and kgutabs to the list
      opt.params <- c(opt.params,
                      c("Fgutabs" = 0.5,
                        "kgutabs" = 2.19)
      )
    }

    #if httk PK params exist, replace the defaults with these
    if(this.dtxsid %in% httk::get_cheminfo(info="dtxsid")){
      httk_params <- httk::parameterize_1comp(dtxsid = DTXSID,
                                              default.to.human = TRUE,
                                              suppress.messages = TRUE,
                                              species = ifelse(tolower(Species) %in% colnames(httk::physiology.data),
                                                               Species,
                                                               "Human"))[c("kelim",
                                                                           "Vdist",
                                                                           "Fgutabs",
                                                                           "kgutabs")]

      if (model %in% "2compartment") {
        httk_params["V1"] <- httk_params["Vdist"]
      }
      replace_names <- intersect(names(opt.params),
                                 names(httk_params))
      opt.params[replace_names] <- httk_params[replace_names]
    }

    #Go parameter by parameter

    #Fgutabs
      if('po' %in% fitdata$Route &
         'iv' %in% fitdata$Route){
        # Need both oral and iv data to get at Fgutabs
        opt.params <- c(opt.params,
                        these.params["Fgutabs"])
        iv.data <- subset(fitdata, Route == "iv" & !is.na(Value))
        iv.data <- subset(iv.data, Time == min(Time))
        iv.data <- subset(iv.data, Dose == max(Dose))
        oral.data <- subset(fitdata,Route == "po" & !is.na(Value))
        oral.data <- subset(oral.data, Value == max(Value))[1, ]
        opt.params["Fgutabs"] <- oral.data$Value / mean(iv.data$Value) * mean(iv.data$Dose) / mean(oral.data$Dose)
        if (is.na(opt.params[["Fgutabs"]]) | opt.params[["Fgutabs"]] > 0){
          opt.params["Fgutabs"] <- 0.5
        }
      } #end if('po' %in% fitdata$Route & 'iv' %in% fitdata$Route)

    #kelim, Vdist/V1. Ralphatobeta
    if ("iv" %in% fitdata$Route){
      #For one moment, assume 1-compartment model describes IV data
      #Conc/Dose = (1/Vdist)*exp(-kel*t)
      #log(Conc/Dose) = log(1/Vdist) + -kel*t
      #of form y = b + m*x, so linear regression works to fit coeffs
      elim.data <- subset(fitdata, Route == "iv" &
                            !is.na(Value) &
                            Dose > 0)
      #normalize concentration to dose
      elim.data$Value <- elim.data$Value / elim.data$Dose

    }else if ("po" %in% fitdata$Route) {
      #if we don't have IV data but we do have oral data
      elim.data <- subset(fitdata, Route == "po" &
                            !is.na(Value) &
                            Dose > 0)
      elim.data$Value <- elim.data$Value/elim.data$Dose
      #get time of Cmax
      max.time <- as.numeric(min(elim.data[Value == max(Value), "Time"]))
      #take only data after the time of Cmax
      elim.data <- subset(elim.data, Time >= max.time)
    } #end if ("iv" %in% fitdata$Route)/else if ("po" %in% fitdata$Route)

    #linear regression of log conc/dose vs. time
    elim_lm <- stats::lm(log(Value) ~ Time,data=elim.data)
    elim_coeff <- coef(elim_lm)
    #take slope
    kelim_tmp <- max(-elim_coeff[2], 0.0001)
    #if kelim_tmp is not NA, then use it.
    if(is.finite(kelim_tmp)){
      opt.params["kelim"] <- kelim_tmp
    }

    #Vdist/V1
    Vd.data <- subset(elim.data, Time == min(Time))
    Vd_tmp <- 1 / mean(Vd.data$Value)
    if(model %in% "1compartment"){
      opt.params["Vdist"] <- Vd_tmp
    }else if(model %in% "2compartment"){
      opt.params["V1"] <- Vd.tmp
    }

    if(model %in% "2compartment"){
      #estimate alpha if possible
      if (length(unique(elim.data$Time)) > 3) {
        #if we have enough time points to separate early and late:
        elim.data.test <- elim.data
        time_early <- sort(unique(elim.data$Time))[3]
        alpha_lm <- stats::lm(log(Value) ~ Time,
                              data =  subset(elim.data,
                                             Time <= time_early)
        )
        alpha <- -coef(alpha_lm)[2]
      }else{
        alpha <- 2
      }

      opt.params["Ralphatokelim"] <- max(alpha / opt.params[["kelim"]], 2)
    }  #end if if(model %in% "2compartment")
  } #end if(model %in% "flat")/else

  #Add a different sigma value for each reference
  refs <- fitdata[,
                  unique(Reference)]
  # Initialize the standard deviations
  #Set a plausible upper bound for sigma:
  TRYSIGMA <- stats::median(fitdata$Value, na.rm = T)

  opt.params[sapply(refs,
                      function(x) paste('sigma2',
                                        x,
                                        sep = '.'))] <- rep(max(TRYSIGMA / 100, 0.1),
                                                            length(refs))

  return(opt.params)
}
