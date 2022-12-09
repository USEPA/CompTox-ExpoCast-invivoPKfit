plot_fit <- function(DTXSID_in,
                     Species_in,
                     Analysis_Type_in,
                     model_in = "winning", #or "all" or by name
                     DF,
                     pk_fit,
                     n = 101,
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


  if(all(Analysis_Type_in %in% "all")){
    Analysis_Type_in <- c("Joint", "Pooled", "Separate")
  }

  DFsub <- DF[DTXSID %in% DTXSID_in &
                Species %in% Species_in]

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

  #Produce a table with unique "experiments":
  #ombinations of chemical, species, reference, route, dose, and medium (blood or plasma).
  DF2 <- DFsub[, .(Max_Time = max(Time),
                Reference = Reference,
                Route = Route,
                Dose = Dose,
                iv = iv,
                Media = Media),
            by = .(DTXSID, Species)]

#Create a column with combined references, to match joint/pooled analyses.
  DF2[, Reference_pooled:=paste(sort(unique(Reference)),
                                collapse = ", "),
      by = .(DTXSID, Species)]

  #Melt to longer format. The effect will be as though we row-bound two versions
  #of DF together: one with the original single referneces, and one with
  #comma-separated combiend references.
  DF3 <- melt(DF2, measure.vars = c("Reference", "Reference_pooled"),
              variable.name = "Reference_type",
              value.name = "Reference")

  #Keep only the unique values of the Reference  column -- there will be
  #duplicated rows for single-reference analyses.
  DF3 <- unique(DF3[, .(DTXSID, Species, Reference,
                        Route, Dose, iv, Media, Max_Time)])




  #subset the fitted parameter data appropriately
  pksub <- subset(pk_fit, DTXSID %in% DTXSID_in &
                    Species %in% Species_in &
                    Analysis_Type %in% Analysis_Type_in)

  #if only winning model selected
  if(all(model_in %in% "winning")){
    pksub <- pksub[winning %in% TRUE, ]
  }

  #if specific model(s) selected by name
  if(!all(model_in %in% c("winning", "all"))){
    pksub <- pksub[model %in% model_in]
  }

  #  Keep only the unique analyses
  pk_wide <- unique(pksub[,
                           .(DTXSID,
                             Species,
                             Analysis_Type,
                             Reference,
                             model,
                             winning)])
#Merge `DF3` and `pk_wide`
  pred_DT <- pk_wide[DF3,
                     on = c("DTXSID", "Species", "Reference"),
                     allow.cartesian = TRUE]

  #Now add time points for prediction.
  pred_DT2 <- pred_DT[, .(Time = 10^(seq(from = log10(0.5/60),
                                         to = log10(Max_Time),
                                         length.out = n))),
                      by = .(DTXSID, Species, Analysis_Type, Reference,
                             model, winning, Route, Dose, iv, Media, Max_Time)]

  #Now evaluate modelfor each analysis & model
  pred_DT2[, Conc := {
    #which model function to use
    if(unique(model) %in% "1compartment"){
      modfun <- "cp_1comp"
    }else if(unique(model) %in% "2compartment"){
      modfun <- "cp_2comp"
    }else{
      modfun <- "cp_flat"
    }
    this_group <- unique(.SD)
    pksub_tmp <- pk_fit[this_group, on = c("DTXSID", "Species", "Analysis_Type", "Reference", "model")]
    #model params
    par <- as.list(pksub_tmp[, `Fitted mean`])
    names(par) <- pksub_tmp[, param_name]

    #call model function
    do.call(modfun,
            args = list(params = par,
                        time = Time,
                        dose = Dose,
                        iv.dose = iv))
  },
  .SDcols = c("DTXSID", "Species", "Analysis_Type", "Reference", "model"),
  by = .(DTXSID, Species, Analysis_Type, Reference, model)]

  #get predicted conc normalized by dose
  pred_DT2[, ConcDose:=Conc/Dose]


  #create a categorical variable for dose
  DFsub[, Dose_cat:=factor(Dose)]

#rename pred_DT2 to predsub
  predsub <- pred_DT2

  #generate plot title
  plot_title <- paste0("DTXSID = ", DTXSID_in, "\n",
                       "Species = ", Species_in, "\n",
                       "Analysis Type =", paste(Analysis_Type_in, collapse= ", "))
  #if plotting a joint analysis (only one winning model for this DTXSID and species),
  #add winning model to the plot title
  if(any(Analysis_Type_in %in% "Joint")){
    plot_title <- paste0(plot_title, "\n",
                         "Winning Model = ", predsub[winning %in% TRUE, unique(model)])
  }

  #plot
  p <- ggplot(data = DFsub,
         aes(x = Time,
             y = ConcDose)) +
    #first plot points with no fill, only stroke
    geom_point(aes(shape = Reference,
                   color = Dose),
               size = 4,
               stroke = 1.5) +
    #then plot points with fill, but alpha mapped to Detect
    geom_point(aes(shape = Reference,
                   color = Dose,
                   fill = Dose,
                   alpha = Detect),
               size = 4) +
    #plot lines for model predictions
    geom_line(data = predsub,
              aes(linetype = model,
                  group = interaction(Analysis_Type, Reference, Dose, model))
    ) +
    facet_grid(rows = vars(Route),
               cols = vars(Media)) +
    scale_color_viridis_c() +
    scale_fill_viridis_c(na.value = NA) +
    scale_shape_manual(values = 21:25) + #use only the 5 shapes where a fill can be added
    #if detect =FALSE, fully transparent; if detect = TRUE, fully solid
    scale_alpha_manual(values = c("Detect" = 1, "Non-Detect" = 0),
                       drop = FALSE) +
    scale_y_log10() +
    annotation_logticks(sides ="l") +
    guides(alpha = guide_legend(override.aes = list(shape = 21,
                                                    color = "black",
                                                    stroke = 1,
                                                    fill = c("black", NA),
                                                    alpha = 1)
    )
    ) +
    ggtitle(plot_title) +
    xlab("Time, hr") +
    ylab("Concentration/Dose") +
    theme_bw()

  #limit y axis scaling to go only 2x smaller than smallest observation
  new_y_min <- DFsub[, min(ConcDose,
                            na.rm = TRUE)/2]
  if(is.finite(new_y_min)){
    p <- p +
      coord_cartesian(ylim = c(new_y_min, NA))
  }

  #log-scale time if so specified
  if(log10_scale_time %in% TRUE){

    p <- p +
      scale_x_log10()  +
      annotation_logticks(sides = "bl")
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
