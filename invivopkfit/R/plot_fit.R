plot_fit <- function(DTXSID_in,
                     Species_in,
                     Analysis_Type_in,
                     model_in = "winning", #or "all" or by name
                     DFsub,
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

  #Get concentration or LOQ
  DFsub[, Conc:=pmax(Value, LOQ, na.rm = TRUE)]
  #Get detection flag
  DFsub[, Detect:=factor(
    ifelse(!is.na(Value), "Detect", "Non-Detect"),
    levels = c("Detect", "Non-Detect")
    )]

  #calculate concentration normalized to dose
  #this allows easier visualization -- all doses plotted together
  DFsub[, ConcDose:=Conc/Dose]

  DFsub[, ConcDose_upper := (Value + Value_SD)/Dose]
  DFsub[, ConcDose_lower := (Value - Value_SD)/Dose]


  #Now add time points for prediction.
  #Interpolate n_interp_time points between each existing time point.
  pred_DT2 <- DFsub[, .(
    Time = {
      timepoints <- unique(Time)
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
                 "Conc", "Detect",
                 "ConcDose", "ConcDose_upper", "ConcDose_lower"))]



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
  pred_DT3[, ConcDose:=Conc/Dose]


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
  p <- ggplot(data = DFsub,
         aes(x = Time,
             y = ConcDose)) +
    geom_blank()

  #add errorbars if so specified
  if(plot_errbars %in% TRUE){
    p <- p +
      geom_errorbar(aes(ymin = ConcDose_lower,
                        ymax = ConcDose_upper,
                        color = Dose))
  }

  #now plot the rest
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
                  group = interaction(Analysis_Type, Studies.Analyzed, Dose, model))
    ) +
    facet_grid(rows = vars(Route),
               cols = vars(Media),
               scales = "free_y") +
    scale_color_viridis_c(name = "Dose, mg/kg") +
    scale_fill_viridis_c(na.value = NA, name = "Dose, mg/kg") +
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
  new_y_min <- DFsub[, min(pmin(ConcDose,
                                ConcDose_lower,
                                na.rm = TRUE),
                            na.rm = TRUE)/2]
  if(!is.finite(new_y_min)) new_y_min <- NA

  new_y_max <- DFsub[, max(pmax(ConcDose,
                                ConcDose_upper,
                                na.rm = TRUE),
                           na.rm = TRUE)*1.05]
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
