

#' Initialize a `pkdata` object
#'
#' Initialize a `pkdata` object
#'
#' A `pkdata` object is a list with elements `data` and `data_info`. `data` is a
#' `data.frame` consisting of a set of concentration vs. time observations to
#' which a single curve will be fitted. `data_info` is a list that provides
#' summary information about the data: DTXSID, Species, routes analyzed, media
#' analyzed, references analyzed, studies analyzed (if "Study" is defined
#' differently from "Reference").
#'
#' The input `data.frame` should contain data for only one `DTXSID` and
#' `Species`. It may contain data for multiple `Route` and/or `Media` values,
#' which can be fitted simultaneously. It may contain data from multiple
#' `Reference` IDs (in the case where a joint or pooled analysis is desired), or
#' from only one `Reference` ID (in the case where a separate analysis is
#' desired for one reference at a time).
#'
#' @param data A `data.frame` -- see Details for required variables. The default
#'   is an empty data frame.
#' @return An object of class `pkdata`

pkdata <- function(data = data.frame(DTXSID = character(),
                                        Compound = character(),
                                        CAS = character(),
                                        Species = character(),
                                        Reference = character(),
                                        Study = character(),
                                        Subject = integer(),
                                        N_Subjects = integer(),
                                        Route = character(),
                                     Dose = numeric(),
                                     Time = numeric(),
                                     Media = character(),
                                     Value = numeric(),
                                     Value_SD = numeric(),
                                     LOQ = numeric()
                                    )
){
  # construct the pkdata object

  #convert data to a data frame
  dat <- as.data.frame(data)

  dat <- initialize_pkdata(dat)

  #get the data info
  #unique DTXSID, Species, References, Studies, Routes
  dat_info <- as.list(unique(dat[c("DTXSID",
                                          "Species")]))

  dat_info$References_Analyzed <- sort(unique(dat$Reference))
  #get a list of studies analyzed
  dat_info$Studies_Analyzed <- sort(unique(dat$Study))
  #get a list of routes analyzed
  dat_info$Routes_Analyzed <- sort(unique(dat$Route))
  #get the number of detects and non-detects by route and medium
  dat_info$n_dat <- aggregate(x = dat$Detect,
                              by = dat[c("Route", "Media")],
                              FUN = function(Detect){
                                c("Detect" = sum(Detect %in% TRUE),
                                      "NonDetect" = sum(Detect %in% FALSE))
                              })
  names(dat_info$n_dat) <- gsub(x = names(dat_info$n_dat),
                                pattern = "Detect.",
                                replacement = "",
                                fixed = TRUE)

  #get empirical tmax
  if(any(dat$Route %in% "po")){
    oral_data <- subset(dat, Route %in% "po")

    oral_peak <- get_peak(x = oral_data$Time,
                          y = log(oral_data[["Conc_Dose"]]))

    dat_info$tmax_oral <- oral_peak$x


  }else{
    dat_info$tmax_oral <- NA_real_
  }

  #absorption and elimination phase points
  #get the number of detects & nondetects before empirical tmax
  dat_info$n_phase <- aggregate(x = dat$Detect,
                              by = list("Phase" = factor(dat$Time <= dat_info$tmax_oral,
                                          levels = c(TRUE, FALSE),
                                          labels = c("absorption", "elimination"))),
                              FUN = function(Detect){
                                c("Detect" = sum(Detect %in% TRUE),
                                  "NonDetect" = sum(Detect %in% FALSE))
                              }
  )
  names(dat_info$n_phase) <- gsub(x = names(dat_info$n_phase),
                                pattern = "Detect.",
                                replacement = "",
                                fixed = TRUE)

  #get time of last detected observation
  dat_info$last_detect_time <- max(dat[dat$Detect %in% TRUE, "Time"])

  #get time of last observation
  dat_info$last_time <- max(dat$Time)

  # construct the pkdata object
  obj <- list("data" = dat,
              "data_info" = dat_info,
  )

  class(obj) <- c(class(obj), "pkdata")

  return(obj)

}

#' Coerce an object to class `pkdata`
#'
#' Given any object that can be coerced to a `data.frame`, convert it to an object of class [pkdata()].
#'
#' @param x The object to coerce to class `pkdata`
#' @return An object of class `pkdata`. This is a `list` with named elements
#'   `data`, `data_info`.
#' @export
#' @author Caroline Ring

as.pkdata.default <- function(x){
  pkdata(as.data.frame(x))
}

#' Check whether an object is of class `pkdata`
#'
#' @param obj The object whose class is to be tested
#' @return TRUE if the object conforms to expectations for class
#'   `pkdata`, FALSE if it does not
#' @export
#' @author Caroline Ring
is.pkdata <- function(obj){
  expected_names <- c("data",
                      "data_info")
  expected_classes <- c("data" = "data.frame",
                        "data_info" = "list")

  check <- inherits(obj, "pkdata") &
    all(expected_names %in% names(obj))

  #if we're OK so far, check whether each named element inherits from the
  #expected class
  if(check %in% TRUE){
    check <- check &
      all(
        sapply(expected_names,
               function(x){
                 inherits(obj[[x]],
                          expected_classes[x])
               })
      )
  }

  return(check)
}

#'Non-compartmental analysis
#'
#'Do non-compartmental analysis of a `pkdata` object.
#'
#'This function calls [PK::nca()] to calculate the following quantities:
#'
#'For intravenous (IV) bolus administration:
#'
#' * `nca.iv.AUC_tlast`: AUC (area under concentration-time curve) evaluated at the last reported time point
#' * `nca.iv.AUC_inf`: AUC evaluated at time t = \eqn{\infty}
#' * `nca.iv.AUMC_inf`: AUMC (area under the first moment curve, i.e. under the AUC vs. time curve) evaluated at time \eqn{t = \infty}
#' * `nca.iv.MRT`: MRT (mean residence time)
#' * `nca.iv.halflife`: Half-life
#' * `nca.iv.Clearance`: Clearance
#' * Vss: apparent volume of distribution at steady-state.
#'
#'For oral bolus administration:
#'
#' * `nca.oral.AUC_tlast`: AUC (area under concentration-time curve) evaluated at the last reported time point
#' * `nca.oral.AUC_inf`: AUC evaluated at time \eqn{t = \infty}
#' * `nca.oral.AUMC_inf`: AUMC (area under the first moment curve, i.e. under the AUC vs. time curve) evaluated at time \eqn{t = \infty}
#' * `nca.oral.MTT`: MTT (mean transit time), the sum of MRT and mean absorption time (MAT)
#'
#'If `obj$data` only contains data for one of those routes, the NCA quantities
#'for the "missing" route will be returned, but filled with `NA_real_.`
#'
#'Additionally, for oral bolus administration, the following quantities are
#'estimated using [get_peak()]:
#'
#' * `nca.oral.tmax`: the time at which peak concentration occurs
#' * `nca.oral.Cmax`: the peak concentration
#'
#'If `dose_norm == TRUE`, then all data are dose-normalized before computing NCA
#'quantities. This means that all dose-dependent NCA quantities reflect an
#'estimate for a 1 mg/kg dose.
#'
#' If `dose_norm == FALSE`, then all quantities are calculated separately for each dose.
#'
#' @param obj A `pkdata` object containing concentration-time data in the element `obj$data`.
#' @param dose_norm Whether to dose-normalize data before performing NCA. Default TRUE.
#' @return A `data.frame` with variables as listed in Details.
nca.pkdata <- function(obj,
                      newdata = NULL,
                      dose_norm = TRUE){

  if(is.null(newdata)){
    newdata <- obj$data
  }

  nca_names <- c("AUC_tlast",
                 "AUC_inf",
                 "AUMC_inf",
                 "MRT",
                 "halflife",
                 "Clearance",
                 "Vss")

  if("iv" %in% newdata$Route){
    iv_data <- subset(newdata,
                      Route %in% "iv")
    if(dose_norm %in% TRUE){
      nca_iv <- do_nca(obs = iv_data,
                       dose_norm = dose_norm)
    }else{
      #do NCA separately for each dose
      dose_list <- split(iv_data,
                         iv_data$Dose)
      nca_iv <- do.call(rbind,
                        sapply(dose_list,
                               function(this_data){
                                 this_nca <- do_nca(obs = this_data,
                                                    dose_norm = FALSE)
                                 this_nca <- cbind(Dose = this_data$Dose)
                                 this_nca
                               })
      )
    }
    names(nca_iv) <- paste0("nca.iv.",
                            names(nca_iv))
  }else{ #if no IV data, fill with NAs
    nca_iv <- data.frame(as.list(rep(NA_real_,
                                     length(nca_names))))
    names(nca_iv) <- paste0("nca.iv.",
                            nca_names)
  }

  if("po" %in% newdata$Route){
    oral_data <- subset(newdata,
                        Route %in% "po")
    if(dose_norm %in% TRUE){
      nca_oral <- get_nca(obs = oral_data,
                          dose_norm = dose_norm)
      nca_oral <- cbind("nca.oral.Dose" = 1,
                        nca_oral)
      max_df <- as.data.frame(
        get_peak(x = oral_data$Time,
                 y = oral_data$Conc/oral_data$Dose)
      )
      names(max_df) <- c("tmax",
                         "Cmax")
      max_df <- cbind("Dose" = 1,
                      max_df)
    }else{
      #do NCA separately for each dose
      dose_list <- split(oral_data,
                         oral_data$Dose)
      nca_oral <- do.call(rbind,
                          sapply(dose_list,
                                 function(this_data){
                                   this_nca <- do_nca(obs = this_data,
                                                      dose_norm = FALSE)
                                   this_nca <- cbind(Dose = this_data$Dose)
                                   this_nca
                                 })
      )
      max_df <- do.call(rbind,
                        sapply(dose_list,
                               function(this_data){
                                 max_list <- as.data.frame(
                                   get_peak(x = this_data$Time,
                                            y = this_data$Conc/this_data$Dose)
                                 )
                                 names(max_list) <- c("tmax",
                                                      "Cmax")
                                 max_list <- cbind("Dose" = this_data$Dose,
                                                   max_list)
                                 max_list
                               }
                        )
      )
    }
    names(nca_oral) <- paste0("nca.oral.",
                              names(nca_oral))

    nca_oral <- merge(nca_oral,
                      max_list,
                      by = "Dose")
  }else{
    nca_oral <- data.frame(as.list(rep(NA_real_,
                                       length(nca_names) + 2)))
    names(nca_oral) <- paste0("nca.oral.",
                              c("Dose",
                                nca_names,
                                c("tmax",
                                  "Cmax")
                              )
    )
  }

  names(nca_oral)[match(c("nca.oral.Clearance",
                          "nca.oral.MRT"))] <- c("nca.oral.Clearance_Fgutabs",
                                                 "nca.oral.MTT")
  #remove oral halflife and Vss estimates because they are not valid
  nca_oral[c("nca.oral.halflife",
             "nca.oral.Vss")] <- NULL

  nca_out <- cbind(nca_iv,
                   nca_oral)
  return(nca_out)
}

#' Plot concentration vs. time data.
#'
#' @param obj A [pkdata()] object with concentration-time data in element `obj$data`.
#' @param newdata Optional: A new set of concentration vs. time data to plot,
#'   different from `obj$data`. Default `NULL` to plot the data in `obj$data`.
#' @return A [ggplot2::ggplot()] plot object.
#' @export
#' @author Caroline Ring
plot_data.pkdata <- function(obj,
                            newdata = NULL,
                            plot_dose_norm = TRUE,
                            plot_log10_conc = FALSE){
  if(is.null(newdata)){
    newdata <- obj$data
  }

  #set the Detect column to a factor with levels "Detect" and "Non-Detect"
  if(!("Detect" %in% names(newdata))){
    newdata$Detect <- TRUE
  }

  newdata$Detect <- factor(newdata$Detect,
                           levels = c(TRUE, FALSE),
                           labels = c("Detect", "Non-Detect"))

  if(!("Conc_SD" %in% names(newdata))){
    newdata$Conc_SD <- 0
  }

  #create generic named columsn for plotting, containing either dose-normalized or
  #non-dose-normalized concentrations, depending on `plot_dose_norm`
  if(plot_dose_norm %in% TRUE){
    newdata$conc_plot <- newdata$Conc_Dose
    newdata$conc_sd_plot <- newdata$Conc_SD_Dose
  }else{
    newdata$conc_plot <- newdata$Conc
    newdata$conc_sd_plot <- newdata$Conc_SD
  }


  #generate plot title
  plot_title <- paste0(obj$data_info$DTXSID,
                       " (", unique(obj$data$Compound), ")\n",
                       "Species = ", obj$data_info$Species, ", ",
                       "Doses = ", paste(signif(
                         sort(unique(obj$data$Dose)),
                         3),
                         collapse = ", "), " mg/kg\n")
  #now plot
  p <- ggplot(newdata,
              aes(x = Time,
                  y = conc_plot)) +
    geom_blank()

  if(plot_log10_conc %in% FALSE){ #plot error bars
    if(plot_dose_norm %in% TRUE){  #map colors to dose
      p <- p +
        geom_errorbar(aes(ymin = conc_plot - conc_sd_plot,
                          ymax = conc_plot + conc_sd_plot,
                          color = Dose))
    }else{   #do not map color to dose
      p <- p +
        geom_errorbar(aes(ymin = conc_plot - conc_sd_plot,
                          ymax = conc_plot + conc_sd_plot))
    }
  } #end if(plot_log10_conc %in% FALSE)

  if(plot_dose_norm %in% TRUE){ #mapping color to dose
    #concentration-dose observation points:
    #shape mapped to Reference, color to Dose, fill yes/no to Detect
    #first plot detected points with white fill
    p <- p +
      #plot detected points with white fill
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference,
                     color = Dose),
                 fill = "white",
                 size = 4,
                 stroke = 1.5) +
      #plot detected points with fill mapped to dose, but alpha mapped to Detect
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference,
                     color = Dose,
                     fill = Dose,
                     alpha = Detect),
                 size = 4,
                 stroke = 1.5) +
      #then plot non-detect points with white fill,
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference,
                      color = Dose),
                  fill = "white",
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect,
                  height = 0) +
      #then plot non-detect points with fill, but alpha mapped to Detect
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference,
                      color = Dose,
                      fill = Dose,
                      alpha = Detect),
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect, height = 0)

    p <- p +
      facet_grid(rows = vars(Route),
                 cols = vars(Media),
                 scales = "free_y")

    p <- p +
      scale_color_viridis_c(name = "Dose, mg/kg") +
      scale_fill_viridis_c(na.value = NA, name = "Dose, mg/kg")

  }else{ #if plot_dose_norm %in% FALSE, don't map color to dose
    #concentration-dose observation points:
    #shape mapped to Reference, fill yes/no to Detect
    #first plot detected points with white fill
    p <- p +
      #plot detected points with white fill
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference),
                 color = "gray50",
                 fill = "white",
                 size = 4,
                 stroke = 1.5) +
      #plot detected points with alpha mapped to Detect
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference,
                     alpha = Detect),
                 color = "gray50",
                 fill = "gray50",
                 size = 4,
                 stroke = 1.5) +
      #then plot non-detect points with white fill,
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference),
                  color = "gray50",
                  fill = "white",
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect,
                  height = 0) +
      #then plot non-detect points with fill, but alpha mapped to Detect
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference,
                      alpha = Detect),
                  color = "gray50",
                  fill = "gray50",
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect, height = 0)

    #do facet_wrap by combinations of route, media, dose
    p <- p + facet_wrap(vars(Route, Media, Dose),
                        labeller = "label_both",
                        scales = "free")
  } #end check if plot_dose_norm %in% TRUE/FALSE

  p <- p +
    scale_shape_manual(values = 21:25) + #use only the 5 shapes where a fill can be added
    #this limits us to visualizing only 5 References
    #but admittedly it's hard to distinguish more than 5 shapes anyway
    #if detect =FALSE, fully transparent; if detect = TRUE, fully solid
    scale_alpha_manual(values = c("Detect" = 1,
                                  "Non-Detect" = 0),
                       drop = FALSE,
                       name = NULL) +
    guides(alpha = guide_legend(override.aes = list(shape = 21,
                                                    color = "black",
                                                    stroke = 1,
                                                    fill = c("black", NA),
                                                    alpha = 1)
    )
    ) +
    labs(title = plot_title) +
    xlab("Time, hr")

  if(plot_dose_norm %in% TRUE){
    p <- p + ylab("Concentration/Dose, (mg/kg)/(mg/L)")
  }else{
    p <- p + ylab("Concentration, mg/L")
  }

  p <- p  +
    theme_bw() +
    theme(plot.title = element_text(size = 12),
          strip.background = element_blank())
  #log-scale y axis if so specified
  if(plot_log10_conc %in% TRUE){
    p <- p + scale_y_log10() + annotation_logticks(sides = "l")
  }

  return(p)

}

initialize_pkdata <- function(data){
  #add any missing columns

  #default: required columns
  dat_default <- data.frame(DTXSID = NA_character_,
                            Compound = NA_character_,
                            CAS = NA_character_,
                            Species = NA_character_,
                            Reference = NA_character_,
                            Study = NA_character_,
                            Subject = NA_integer_,
                            N_Subjects = NA_integer_,
                            Route = NA_character_,
                            Dose = NA_real_,
                            Time = NA_real_,
                            Media = NA_character_,
                            Value = NA_real_,
                            Value_SD = NA_real_,
                            LOQ = NA_real_)

  missing_cols <- setdiff(names(dat_default),
                          names(dat))

  if(length(missing_cols)>0){
    dat <- cbind(dat,
                 dat_default[missing_cols])
  }


  #create detect column if not there already
  if(!("Detect" %in% names(dat))){
    dat$Detect <- !is.na(dat$Value)
  }

  #create Conc column if not there already
  if(!("Conc" %in% names(dat))){
    dat$Conc <- pmax(dat$Value, dat$LOQ, na.rm = TRUE)
  }

  #get SDs if not there already
  if(!("Conc_SD" %in% names(dat)) &
     "Value_SD" %in% names(dat)){
    dat$Conc_SD <- dat$Value_SD

  }

  #create dose-normalized variables
  dat$Conc_Dose <- dat$Conc/dat$Dose
  dat$Value_Dose <- dat$Value/dat$Dose
  dat$Value_SD_Dose <- dat$Value_SD/dat$Dose
  dat$LOQ_Dose <- dat$LOQ/dat$Dose
  dat$Conc_SD_Dose <- dat$Conc_SD/dat$Dose
}
