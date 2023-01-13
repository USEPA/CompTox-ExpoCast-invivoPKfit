#'Plot TK model fits
#'
#'Plot fitted toxicokinetic model curves with concentration vs. time data
#'
#'This function can plot dose-normalized data, non-dose-normalized data, can
#'split doses into separate facets or not.
#'@param DTXSID_in DSSTox Substance ID of chemical to be plotted.
#'@param Species_in Species in data to be plotted
#'@param Analysis_Type_in "Joint", "Pooled", or "Separate" (whether different
#'  studies were fit together with separate error SDs; together with one error
#'  SD; or separately)
#'@param model_in Which model or models to plot fits for. Options include any or
#'  all of "1compartment", "2compartment", and "flat"; "winning" to plot only
#'  the winning model (whichever had lowest AIC); or "all" as a short form to
#'  plot all three models.
#'@param DFsub Pre-processed concentration vs. time data merged with fitted PK
#'  parameters. (I will create a separate function to produce this)
#'@param plot_dose_norm Logical: TRUE to plot dose-normalized concentration vs.
#'  time data and model fits. FALSE to plot non-dose-normalized data and fits.
#'  Default TRUE. If `plot_dose_norm == FALSE`, it's typically recommended to
#'  use `split_dose == TRUE`, because non-dose-normalized concentrations for
#'  different dose groups will likely be on different scales.
#'@param split_dose Logical: TRUE to split separate dose groups into separate
#'  plot facets. FALSE to facet only by Route and Media. If `split_dose ==
#'  FALSE`, then dose groups will be distinguished using color.
#'@param n_interp_time The number of time points to interpolate between each
#'  existing time point in the data set for plotting model fits. Default is 10.
#'@param limit_y_axis Whether to limit the scale of the y axis to only a small
#'  distance above and below the min/max of the concentration vs. time data.
#'  Default FALSE. (You may want to try TRUE if you are log-scaling the y axis.)
#'@param log10_scale_y Logical: Whether to log-scale the y axis. Default FALSE.
#'@param plot_errbars Logical: Whether to plot sample errorbars for
#'  concentration vs. time data (when the observations have associated sample
#'  standard deviations). Default TRUE. You may want to try FALSE if you are
#'  log-scaling the y axis.
#'@param log10_scale_time Logical: Whether to log-scale the x-axis (time).
#'  Default FALSE.
#'@param save_plot Logical: Whether to save the plot as a file. Default TRUE.
#'@param return_plot Logical: Whether to return the plot object. Default FALSE.
#'@param file_path Path to save file if `save_plot == TRUE`. Default
#'  `"inst/ext/plots"`.
#'@param file_suffix Suffix to append to the file name, if any. Default NULL for
#'  no suffix.
#'@param file_format File format to save the plot. Default `"pdf"`.
#'@param ggsave_args A named list of additional arguments to `ggsave()` to use
#'  when saving the plot. Default
#'   ```
#'   list(scale = 1,
#'width = 9,
#'height = 10,
#'units = "in",
#'dpi = 300,
#'limitsize = TRUE,
#'bg = NULL)
#'```
#'@param verbose Logical: Whether to print verbose messages about plotting
#'  progress. Default TRUE.
#'@return If `return_plot == TRUE`, a `ggplot2` object. If `return_plot ==
#'  FALSE`, returns 0.
#' @author Caroline Ring
#' @export
#'
plot_fit <- function(DTXSID_in,
                     Species_in,
                     Analysis_Type_in,
                     model_in = "winning", #or "all" or by name
                     DFsub,
                     plot_dose_norm = TRUE,
                     split_dose = FALSE,
                     n_interp_time = 10,
                     limit_y_axis = FALSE,
                     log10_scale_y = FALSE,
                     plot_errbars = TRUE,
                     log10_scale_time = FALSE,
                     save_plot = TRUE,
                     return_plot = FALSE,
                     file_path = "inst/ext/plots/",
                     file_suffix = NULL,
                     file_format = "pdf",
                     ggsave_args = list(scale = 1,
                                        width = 9,
                                        height = 10,
                                        units = "in",
                                        dpi = 300,
                                        limitsize = TRUE,
                                        bg = NULL),

                     verbose = TRUE){

  if(verbose %in% TRUE){
    message(paste0("Plotting:\n",
      "DTXSID = ", DTXSID_in, "\n",
                  "Species = ", Species_in, "\n",
                  "Analysis Type = ", Analysis_Type_in, "\n",
                  "models = ", paste(model_in, collapse = ", ")))
  }


  if(all(Analysis_Type_in %in% "all")){
    Analysis_Type_in <- c("Joint", "Pooled", "Separate")
  }

  if(any(model_in %in% "all")){
    model_in <- c("flat", "1compartment", "2compartment")
  }

  DFsub <- DFsub[Analysis_Type %in% Analysis_Type_in, ]

 #calculate lower and upper bounds for obs with sample SD
  DFsub[, Conc_Dose_upper := (Value + Value_SD)/Dose]
  DFsub[, Conc_Dose_lower := (Value - Value_SD)/Dose]

  DFsub[, Conc_upper := (Value + Value_SD)]
  DFsub[, Conc_lower := (Value - Value_SD)]

  #ensure that variable Detect is a factor with levels "Detect" and "Non-Detect"
  DFsub[, Detect:=factor(Detect,
                         levels = c("Detect",
                                    "Non-Detect"))]

  #Now add time points for prediction.
  #Interpolate n_interp_time points between each existing time point.
  pred_DT2 <- DFsub[, .(
    Time = {
      timepoints <- sort(unique(c(0,Time)))
      time_interp <- sapply(1:(length(timepoints)-1),
                            function(i){
                              seq(from = timepoints[i],
                                  to = timepoints[i+1],
                                  length.out = n_interp_time)
                            }
      )
      sort(unique(time_interp))
    }
  ),
  by = setdiff(names(DFsub),
               c("Time", "Value", "Value_SD", "LOQ", "N_Subjects",
                 "Conc", "Detect", "Conc_upper", "Conc_lower",
                 "Conc_Dose", "Conc_Dose_upper", "Conc_Dose_lower",
                 "Value_Dose", "Value_SD_Dose", "LOQ_Dose"))]



  #Now evaluate model function for each analysis & model
  pred_DT2[, Conc.flat := {
    params <- as.list(unique(.SD))
    names(params) <- gsub(names(params),
                          pattern = ".flat",
                          replacement = "",
                          fixed = TRUE)
    #drop any NA params
    params <- params[sapply(params,
                            is.finite)]
    #if all params are missing, return NA;
    #otehrwise, evaluate model
    if(length(params)>0){
    if(!("Rblood2plasma" %in% names(params))){
      params$Rblood2plasma <- 1
    }
      cp_flat(time = Time,
               params = params,
               dose = Dose,
               iv.dose = iv,
               medium = Media)
    }else{
      NA_real_
    }
    },
           .SDcols= c("Vdist.flat",
                      "Fgutabs.flat",
                      "Fgutabs_Vdist.flat",
                      "Rblood2plasma.flat"),
           by = .(Analysis_Type,
                  DTXSID,
                  Species,
                  References.Analyzed,
                  Studies.Analyzed)]

  pred_DT2[, Conc.1compartment := {
    params <- as.list(unique(.SD))
    names(params) <- gsub(names(params),
                          pattern = ".1compartment",
                          replacement = "",
                          fixed = TRUE)

    #drop any NA params
    params <- params[sapply(params,
                            is.finite)]
    #if all params are missing, return NA;
    #otehrwise, evaluate model
    if(length(params)>0){
      if(!("Rblood2plasma" %in% names(params))){
        params$Rblood2plasma <- 1
      }
      cp_1comp(time = Time,
               params = params,
               dose = Dose,
               iv.dose = iv,
               medium = Media)
    }else{
      NA_real_
    }
  },
  .SDcols= c("Vdist.1compartment",
             "Fgutabs.1compartment",
             "Fgutabs_Vdist.1compartment",
             "kelim.1compartment",
             "kgutabs.1compartment",
             "Rblood2plasma.1compartment"),
  by = .(Analysis_Type,
         DTXSID,
         Species,
         References.Analyzed,
         Studies.Analyzed)]

  pred_DT2[, Conc.2compartment := {
    params <- as.list(unique(.SD))
    names(params) <- gsub(names(params),
                          pattern = ".2compartment",
                          replacement = "",
                          fixed = TRUE)

    #drop any NA params
    params <- params[sapply(params,
                            is.finite)]
    #if all params are missing, return NA;
    #otehrwise, evaluate model
    if(length(params)>0){
      if(!("Rblood2plasma" %in% names(params))){
        params$Rblood2plasma <- 1
      }
    cp_2comp(time = Time,
                    params = params,
                    dose = Dose,
                    iv.dose = iv,
                    medium = Media)
    }else{
      NA_real_
    }
  },
  .SDcols= c("V1.2compartment",
             "Fgutabs.2compartment",
             "Fgutabs_V1.2compartment",
             "kelim.2compartment",
             "kgutabs.2compartment",
             "k12.2compartment",
             "k21.2compartment",
             "Rblood2plasma.2compartment"),
  by = .(Analysis_Type,
         DTXSID,
         Species,
         References.Analyzed,
         Studies.Analyzed)]

  #melt pred_DT2
  pred_DT3 <- melt(pred_DT2,
                   id.vars = c(idcols, "winmodel", "Time",
                               "Route", "Dose", "iv", "Media"),
                   measure.vars = c("Conc.flat",
                                    "Conc.1compartment",
                                    "Conc.2compartment"),
                   variable.name = "model",
                   value.name = "Conc")
  pred_DT3[, model := gsub(x = model,
                           pattern= "Conc.",
                           replacement = "",
                           fixed = TRUE)]

  pred_DT3[model == winmodel, winning := TRUE]
  pred_DT3[model != winmodel, winning := FALSE]

  #keep only the model(s) specified
  if(all(model_in %in% "winning")){
    pred_DT3 <- pred_DT3[winning == TRUE,]
  }else if(!(any(model_in %in% "winning"))){
    pred_DT3 <- pred_DT3[model %in% model_in]
  }else{ #if model_in is "winning" and something else
    pred_DT3 <-  pred_DT3[winning == TRUE |
                            model %in% model_in,]
  }


  #get predicted conc normalized by dose
  pred_DT3[, Conc_Dose:=Conc/Dose]


  #create a categorical variable for dose
  DFsub[, Dose_cat:=factor(Dose)]

#rename pred_DT2 to predsub
  predsub <- pred_DT3

  #generate plot title
  plot_title <- paste0(DTXSID_in,
                       " (", DFsub[, unique(Compound)], ")\n",
                       "Species = ", Species_in, ", ",
                       "Doses = ", paste(signif(
                         sort(unique(DFsub$Dose)),
                         3),
                                         collapse = ", "), " mg/kg\n",
  "Analysis Type = ", paste(Analysis_Type_in, collapse= ", "))

  #if plotting a joint analysis (only one winning model for this DTXSID and species),
  #add winning model to the plot title
  if(any(Analysis_Type_in %in% "Joint")){
    winning_model <- predsub[winning %in% TRUE, unique(model)]
    plot_subtitle <- paste0("Winning Model = ",
                         winning_model)
    #get coeffs of winning model
    par_names <- grep(x = names(DFsub),
                pattern = paste0(".", winning_model),
                fixed = TRUE,
                value = TRUE)
    par <- unique(DFsub[, .SD, .SDcols = par_names])
    setnames(par,
             names(par),
             gsub(x = names(par),
                  pattern = paste0(".", winning_model),
                  replacement = ""))
    par <- unlist(par)
    #remove any NA or infinite parameters
    par <- par[is.finite(par)]

    #paste into a comma-separated list
    par_char <- paste(
      paste(names(par),
            signif(par, 3),
            sep=" = "),
      collapse = ", ")

    plot_subtitle <- paste0(plot_subtitle, "\n",
                         par_char)
  }else{
    plot_subtitle <- NULL
  }

  #plot
  if(plot_dose_norm %in% TRUE){
    p <- ggplot(data = DFsub,
                aes(x = Time,
                    y = Conc_Dose)) +
      geom_blank()

    #add errorbars if so specified
    if(plot_errbars %in% TRUE){
      if(split_dose %in% FALSE){
        p <- p +
          geom_errorbar(aes(ymin = Conc_Dose_lower,
                            ymax = Conc_Dose_upper,
                            color = Dose))
      }else{ #if split_dose == TRUE, don't map color to Dose
        p <- p + geom_errorbar(aes(ymin = Conc_Dose_lower,
                                 ymax = Conc_Dose_upper))
      }
    }
  }else{ #if plot_dose_norm == FALSE
    p <- ggplot(data = DFsub,
                aes(x = Time,
                    y = Conc)) +
      geom_blank()

    #add errorbars if so specified
    if(plot_errbars %in% TRUE){
      if(split_dose %in% FALSE){
        p <- p +
          geom_errorbar(aes(ymin = Conc_lower,
                            ymax = Conc_upper,
                            color = Dose))
      }else{ #if split_dose %in% TRUE, don't map color to dose
        p <- p +
          geom_errorbar(aes(ymin = Conc_lower,
                            ymax = Conc_upper))
      }
    }
  }

  #now plot the rest


    if(split_dose %in% FALSE){
      #concentration-dose observation points:
      #shape mapped to Reference, color to Dose, fill yes/no to Detect
      #first plot points with white fill
      p <- p + geom_point(aes(shape = Reference,
                              color = Dose),
                          fill = "white",
                          size = 4,
                          stroke = 1.5) +
        #then plot points with fill, but alpha mapped to Detect
        geom_point(aes(shape = Reference,
                       color = Dose,
                       fill = Dose,
                       alpha = Detect),
                   size = 4,
                   stroke = 1.5) +
        #plot lines for model predictions
        geom_line(data = predsub,
                  aes(linetype = model,
                      group = interaction(Analysis_Type,
                                          Studies.Analyzed,
                                          Dose,
                                          model))
        )

      p <- p +
    facet_grid(rows = vars(Route),
               cols = vars(Media),
               scales = "free_y")

      p <- p +
        scale_color_viridis_c(name = "Dose, mg/kg") +
        scale_fill_viridis_c(na.value = NA, name = "Dose, mg/kg")
    }else{ #if split_dose == TRUE, don't map color to Dose
      #concentration-dose observation points:
      #shape mapped to Reference, fill yes/no to Detect
      #first plot points with white fill
      p <- p + geom_point(aes(shape = Reference),
                          color = "gray50",
                          fill = "white",
                          size = 4,
                          stroke = 1.5) +
        #then plot points with fill, but alpha mapped to Detect
        geom_point(aes(shape = Reference,
                       alpha = Detect),
                   color = "gray50",
                   fill =  "gray50",
                   size = 4,
                   stroke = 1.5) +
        #plot lines for model predictions
        geom_line(data = predsub,
                  aes(linetype = model,
                      group = interaction(Analysis_Type,
                                          Studies.Analyzed,
                                          Dose,
                                          model))
        )
      p <- p + facet_wrap(vars(Route, Media, Dose),
                          labeller = "label_both",
                          scales = "free_y")
    }

    p <- p +
    scale_shape_manual(values = 21:25) + #use only the 5 shapes where a fill can be added
    #this limits us to visualizing only 5 References
    #but admittedly it's hard to distinguish more than 5 shapes anyway
    #if detect =FALSE, fully transparent; if detect = TRUE, fully solid
    scale_alpha_manual(values = c("Detect" = 1, "Non-Detect" = 0),
                       drop = FALSE,
                       name = NULL) +
    guides(alpha = guide_legend(override.aes = list(shape = 21,
                                                    color = "black",
                                                    stroke = 1,
                                                    fill = c("black", NA),
                                                    alpha = 1)
    )
    ) +
    labs(title = plot_title,
            subtitle = plot_subtitle) +
    xlab("Time, hr") +
    ylab("Concentration/Dose") +
    theme_bw() +
      theme(plot.title = element_text(size = 12),
            strip.background = element_blank())


  #limit y axis scaling to go only 2x smaller than smallest observation
    if(limit_y_axis %in% TRUE){
      if(plot_dose_norm %in% TRUE){
  new_y_min <- DFsub[, min(pmin(Conc_Dose,
                                Conc_Dose_lower,
                                na.rm = TRUE),
                            na.rm = TRUE)/2]

  new_y_max <- DFsub[, max(pmax(Conc_Dose,
                                Conc_Dose_upper,
                                na.rm = TRUE),
                           na.rm = TRUE)*1.05]
      }else{
        new_y_min <- DFsub[, min(pmin(Conc,
                                      Conc_lower,
                                      na.rm = TRUE),
                                 na.rm = TRUE)/2]


        new_y_max <- DFsub[, max(pmax(Conc,
                                      Conc_upper,
                                      na.rm = TRUE),
                                 na.rm = TRUE)*1.05]
      }
      if(!is.finite(new_y_min)) new_y_min <- NA
  if(!is.finite(new_y_max)) new_y_max <- NA

    p <- p +
      coord_cartesian(ylim = c(new_y_min, new_y_max))
    }


  #log-scale concentration/dose if so specified
  if(log10_scale_y %in% TRUE){
    p <- p + scale_y_log10()
  }

  #log-scale time if so specified
  if(log10_scale_time %in% TRUE){

    p <- p +
      scale_x_log10()
  }

  #add logtick annotations if appropriate
  if(log10_scale_y %in% TRUE &
     log10_scale_time %in% TRUE){
    p <- p + annotation_logticks(sides = "bl")
  }else if(log10_scale_y %in% TRUE &
           log10_scale_time %in% FALSE){
    p <- p + annotation_logticks(sides = "l")
  }else if(log10_scale_y %in% FALSE &
           log10_scale_time %in% TRUE){
    p <- p + annotation_logticks(sides = "b")
  }

#set up filename for saving
  filename <- paste0(DTXSID_in,
                     "_",
                     Species_in,
                     "_",
                     paste(Analysis_Type_in, collapse = "_"),
                     "_",
                     paste(model_in, collapse = "_"))

  if(length(file_suffix) > 0){
    filename <- paste0(filename,
                       "_",
                       file_suffix)
  }

  filename <- paste0(filename, ".", file_format)


if(save_plot %in% TRUE){

  #and save
  do.call(ggsave,
          args = c(list(filename = filename,
                        plot = p,
                        device = file_format,
                        path = file_path),
                   ggsave_args)
          )
}

  if(return_plot %in% TRUE){
    return(p)
  }else{
    return(0)
  }

}
