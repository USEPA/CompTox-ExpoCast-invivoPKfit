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
#'@param obs_data A `data.table`: Observed concentration vs. time data to be
#'  plotted, pre-processed as from [preprocess_data()]. Must have variables
#'  `DTXSID`, `Species`, `Route`, `Media`, `iv`, `Dose`, `Time`, `Conc`,
#'  `Detect`, `Conc_SD`.
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
#'@author Caroline Ring
#'@export
#'
plot_fit <- function(fit_flat,
                     fit_1comp,
                     fit_2comp,
                     DTXSID_in,
                     Species_in,
                     Analysis_Type_in,
                     Studies.Analyzed_in,
                     model_in = "winning", #or "all" or by name
                     fit_log_conc,
                     fit_conc_dose,
                     rescale_time,
                     obs_data,
                     plot_dose_norm = TRUE,
                     split_dose = FALSE,
                     n_interp_time = 10,
                     limit_y_axis = FALSE,
                     log10_scale_y = FALSE,
                     plot_errbars = TRUE,
                     log10_scale_time = FALSE,
                     jitter_nondetect = 0.1,
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

  pk_fit <- merge_fits(fit_flat = fit_flat,
                       fit_1comp = fit_1comp,
                       fit_2comp = fit_2comp)

  if(all(Analysis_Type_in %in% "all")){
    Analysis_Type_in <- c("Joint", "Pooled", "Separate")
  }

  if(any(model_in %in% "all")){
    model_in <- c("flat", "1compartment", "2compartment")
  }

  if(verbose %in% TRUE){
    message(paste0("Plotting:\n",
      "DTXSID = ", DTXSID_in, "\n",
                  "Species = ", Species_in, "\n",
                  "Analysis Type = ", paste(Analysis_Type_in, collapse = ", "), "\n",
                  "models = ", paste(model_in, collapse = ", ")))
  }


 #calculate lower and upper bounds for obs with sample SD
  obs_data[, Conc_Dose_upper := (Value + Value_SD)/Dose]
  obs_data[, Conc_Dose_lower := (Value - Value_SD)/Dose]

  obs_data[, Conc_upper := (Value + Value_SD)]
  obs_data[, Conc_lower := (Value - Value_SD)]

  #ensure that variable Detect is a factor with levels "Detect" and "Non-Detect"
  obs_data[, Detect:=factor(Detect,
                         levels = c("Detect",
                                    "Non-Detect"))]

  #Now create interpolated time points for prediction.
  #Interpolate n_interp_time points between each existing time point.
  #Go by DTXSID, Species, Route, Media, Dose, Subject.
  pred_DT2 <- obs_data[, .(
    Time = {
      timepoints <- sort(unique(c(0,0,Time)))
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
  by = .(DTXSID, Compound, Species, Route, iv, Media, Study, Dose)]


  #Now evaluate model function for pred_DT2, for each analysis & model Set up a
  #table of all combinations of DTXSID, Species, analsis type, model, and
  #Studies.Analyzed to be plotted
  pk_sub <- pk_fit[DTXSID %in% DTXSID_in &
                     Species %in% Species_in &
                     Analysis_Type %in% Analysis_Type_in]

  if(all(model_in %in% "winning")){
    tmp <- unique(pk_sub[, .(DTXSID, Species, Analysis_Type, Studies.Analyzed, winmodel)])
    setnames(tmp,
             "winmodel",
             "model")
  }else{
  tmp <- unique(pk_sub[, .(DTXSID, Species, Analysis_Type, Studies.Analyzed)])
  tmp <- tmp[, .(model = model_in), by = names(tmp)]
  }

  setnames(tmp,
           names(tmp),
           paste(names(tmp), "in", sep = "_"))

#Now use data.table syntax to loop over the rows of `tmp` and evaluate predictions for each one
  pred_data <- tmp[,
                  {
                    studies_list <- strsplit(Studies.Analyzed_in,
                                            split = ", ")[[1]]
                    #get data.table of newdata to be predicted
                    #leave out DTXSID, Species columns
                    #since these will be added back in from fit_combs
                    newdata <- pred_DT2[DTXSID %in% DTXSID_in &
                                         Species %in% Species_in &
                                         Study %in% studies_list,
                                       .SD,
                                       .SDcols = setdiff(names(pred_DT2),
                                                         c("DTXSID", "Species"))]
                    get_predictions(pk_fit = pk_fit,
                              newdata = newdata,
                              DTXSID_in = DTXSID_in,
                              Species_in = Species_in,
                              Analysis_Type_in = Analysis_Type_in,
                              Studies.Analyzed_in = Studies.Analyzed_in,
                              model_in = model_in)
                  }, by = .(DTXSID_in,
                          Species_in,
                          Analysis_Type_in,
                          Studies.Analyzed_in,
                          model_in)
  ]

  setnames(pred_data,
           "pred_conc",
           "Conc")

  setnames(pred_data,
           names(pred_data),
           gsub(x = names(pred_data),
                pattern= "_in",
                replacement = ""))

  #get predicted conc normalized by dose
  pred_data[, Conc_Dose:=Conc/Dose]

  #create a categorical variable for dose in observed data
  obs_data[, Dose_cat:=factor(Dose)]

  #generate plot title
  plot_title <- paste0(DTXSID_in,
                       " (", obs_data[, unique(Compound)], ")\n",
                       "Species = ", Species_in, ", ",
                       "Doses = ", paste(signif(
                         sort(unique(obs_data$Dose)),
                         3),
                                         collapse = ", "), " mg/kg\n",
  "Analysis Type = ", paste(Analysis_Type_in, collapse= ", "), "\n",
  "Fitting options: ", "log-transform ", fit_log_conc,
  "; dose-normalize ", fit_conc_dose,
  "; rescale time ", rescale_time)

  #if plotting a joint analysis (only one winning model for this DTXSID and species),
  #add winning model to the plot title
  if(any(Analysis_Type_in %in% "Joint")){
    #get winning model for this dataset & joint analysis
    pk_sub <- pk_fit[DTXSID %in% DTXSID_in &
                       Species %in% Species_in &
                       Analysis_Type %in% "Joint", ]
    winning_model <- pk_sub[, winmodel]
    plot_subtitle <- paste0("Winning Model = ",
                         winning_model)
if(!(grepl(x = winning_model, pattern = "None"))){
    #get coeffs of winning model
    #get parameter names for this model
    param_names <- get_model_paramnames(model = winning_model)
    #get names of appropriate columns of pk_sub -- named [param].[model]
    param_names_dot <- paste(param_names, winning_model, sep = ".")
    #extract the appropriate columns from pk_sub
    par <- unique(pk_sub[, .SD, .SDcols = param_names_dot])
    #remove the ".[model]" suffix
    setnames(par,
             names(par),
             gsub(x = names(par),
                  pattern = paste0(".", winning_model),
                  replacement = ""))
    #convert from data.table to vector
    par <- unlist(par)
    #remove any NA or infinite parameters
    par <- par[is.finite(par)]
    #paste into a comma-separated list
    par_char <- paste(
      paste(names(par),
            signif(par, 3), #keep 3 sigfigs
            sep=" = "),
      collapse = ", ")

    plot_subtitle <- paste0(plot_subtitle, "\n",
                         par_char)
  }else{ #if winning model is "None (no fit)"
    plot_subtitle <- NULL
  }
  }else{ #if analysis type is not "Joint"
    plot_subtitle <- NULL
  }

  #plot
  if(plot_dose_norm %in% TRUE){
    p <- ggplot(data = obs_data,
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
    p <- ggplot(data = obs_data,
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
      #first plot detected points with white fill
      p <- p +
        geom_point(data = obs_data[Detect %in% "Detect"],
                   aes(shape = Reference,
                              color = Dose),
                          fill = "white",
                          size = 4,
                          stroke = 1.5) +
        #then plot detected points with fill, but alpha mapped to Detect
        geom_point(data = obs_data[Detect %in% "Detect"],
                   aes(shape = Reference,
                       color = Dose,
                       fill = Dose,
                       alpha = Detect),
                   size = 4,
                   stroke = 1.5) +
        #then plot non-detect points with white fill,
        #and position jittered
        geom_jitter(data = obs_data[Detect %in% "Non-Detect"],
                    aes(shape = Reference,
                        color = Dose),
                    fill = "white",
                    size = 4,
                    stroke = 1.5,
                    width = jitter_nondetect,
                    height = 0) +
        #then plot non-detect points with fill, but alpha mapped to Detect
        #and position jittered
        geom_jitter(data = obs_data[Detect %in% "Non-Detect"],
                   aes(shape = Reference,
                       color = Dose,
                       fill = Dose,
                       alpha = Detect),
                   size = 4,
                   stroke = 1.5,
                   width = jitter_nondetect, height = 0) +
        #plot lines for model predictions
        geom_line(data = pred_data,
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
        geom_line(data = pred_data,
                  aes(linetype = model,
                      group = interaction(Analysis_Type,
                                          Studies.Analyzed,
                                          Dose,
                                          model))
        )
      p <- p + facet_wrap(vars(Route, Media, Dose),
                          labeller = "label_both",
                          scales = "free")
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


  #limit y axis scaling to go only 2x smaller than smallest observation
    if(limit_y_axis %in% TRUE){
      if(plot_dose_norm %in% TRUE){
  new_y_min <- obs_data[, min(pmin(Conc_Dose,
                                Conc_Dose_lower,
                                na.rm = TRUE),
                            na.rm = TRUE)/2]

  new_y_max <- obs_data[, max(pmax(Conc_Dose,
                                Conc_Dose_upper,
                                na.rm = TRUE),
                           na.rm = TRUE)*1.05]
      }else{
        new_y_min <- obs_data[, min(pmin(Conc,
                                      Conc_lower,
                                      na.rm = TRUE),
                                 na.rm = TRUE)/2]


        new_y_max <- obs_data[, max(pmax(Conc,
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
