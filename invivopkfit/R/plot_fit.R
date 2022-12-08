plot_fit <- function(DTXSID_in,
                     Species_in,
                     Analysis_Type_in,
                     Reference_in,
                     model_in = c("1compartment", "2compartment", "flat"),
                     fitdata,
                     pk_fit,
                     save_plot = TRUE,
                     return_plot = FALSE,
                     file_path = "inst/ext/plots/",
                     file_suffix = NULL,
                     file_format = "pdf",
                     ggsave_args = list(scale = 1,
                                        width = 8.5,
                                        height = 11,
                                        units = "in",
                                        dpi = 300,
                                        limitsize = TRUE,
                                        bg = NULL)){
  #if more than one refrence in comma separated list, get a vector of references
  reflist <- strsplit(Reference_in,
                      split =  ", ")[[1]]

#set up filename for saving
  filename <- paste0(DTXSID_in,
                     "_",
                     Species_in,
                     "_",
                     Reference_in,
                     "_",
                     Analysis_Type_in,
                     "_",
                     paste(model_in, collapse = "_"))

  if(length(file_suffix) > 0){
    filename <- paste0(filename,
                       "_",
                       file_suffix)
  }

  filename <- paste0(filename, ".", file_format)

  # set up plot title
 plot_title <- paste0(DTXSID_in, "\n",
         "Species = ", Species_in, "\n",
         "Reference = ", Reference_in, "\n",
         "Analysis Type = ", Analysis_Type_in,
         "Model = ", paste(model_in, collapse = ", "))

  #subset the data appropriately
  fitsub <- subset(fitdata, DTXSID %in% DTXSID_in &
                     Species %in% Species_in &
                     Reference %in% reflist)

  #subset the fitted parameter data appropriately
  pksub <- subset(pk_fit, DTXSID %in% DTXSID_in &
                    Species %in% Species_in &
                    Reference %in% Reference_in &
                    Analysis_Type %in% Analysis_Type_in)

  #if user specified the winning model, figure out which one that is
  if(all(model_in %in% "winning")){
    #get models & AICs
    pk_aic <- unique(pksub[, .(model, AIC)])
    #sort models by ascending AIC
    setorder(pk_aic, AIC)
    #winning model is the first row
    model_in <- pk_aic[1, AIC]
  }

  if(all(model_in %in% "all")){
    model_in = c("1compartment", "2compartment", "flat")
  }

  #now subset down to the selected model or models
  pksub <- pksub[model %in% model_in]

  #Get concentration or LOQ
  fitsub[, Conc:=pmax(Value, LOQ, na.rm = TRUE)]
  #Get detection flag
  fitsub[, Detect:=!is.na(Value)]

  #calculate concentration normalized to dose
  #this allows easier visualization -- all doses plotted together
  fitsub[, ConcDose:=Conc/Dose]

  #Get a list of lists of model params, one for each model_in
  params <- sapply(model_in,
                   function(mod) {
                     parval <- pksub[model %in% mod,
                                     `Fitted mean`]
                     names(parval) <- pksub[model %in% mod,
                                            param_name]
                     return(as.list(parval))
                   },
                   USE.NAMES = TRUE,
                   simplify = FALSE)

  #color will represent route (po/iv)
  #filled symbols will represent detects (>LOQ)
  #emtpy symbols will represent nondetects (<LOQ)
  #So we need to map shape to ND
  shapevect <- c(
                 "FALSE" = 1, #open circle
                 "TRUE" = 19 #filled circle
  )

  colorvect <- RColorBrewer::brewer.pal(n = 3, name = "Set2")[1:2]
  names(colorvect) <- c("po", "iv")


  #plot the conc-time-dose data
  p <- ggplot(data = fitsub,
              aes(x = Time,
                  y = ConcDose,
                  color = Route,
                  shape = Detect
              ),
              size = 3) +
    geom_point()

  #get a table of unique combinations of model & route
  modelroute <- expand.grid(model = model_in,
                            route = unique(fitsub$Route),
                            stringsAsFactors = FALSE)

  #For each model & route, add a stat_function() layer to plot model predictions
  #Produce a list of stat_function() layers

  layer_list <- mapply(function(mod, route){
    #which model function to use
    if(mod %in% "1compartment"){
      modfun <- "cp_1comp"
    }else if(mod %in% "2compartment"){
      modfun <- "cp_2comp"
    }else{
      modfun <- "cp_flat"
    }

    stat_function(fun = function(x, params, route) do.call(modfun,
                                                        args = list("params" = params,
                                                                    "time" = x,
                                                                    "dose" = 1,
                                                                    "iv" = route %in% "iv")),
                  do.call(aes, args = list(linetype = mod,
                                           color = route)),
                  args = list("params" = params[[mod]],
                              "route" = route),
                  geom = "line")

  },
  mod = modelroute$model,
  route = modelroute$route,
  SIMPLIFY = FALSE)


  #now add all the layers
  p_all <- Reduce(`+`,
                  c(list(p), layer_list))

  #add scales
  p_all <- p_all +
    scale_y_log10() +
    annotation_logticks(sides = "l") +
    scale_shape_manual(values = shapevect,
                       name = NULL,
                       breaks = c("TRUE",
                                  "FALSE"),
                       labels = c("Detect",
                                  "Nondetect"),
                       drop = TRUE) +
    scale_linetype_manual(name = "Model",
                          values = c("1compartment" = 1,
                                     "2compartment" = 2,
                                     "flat" = 3),
                          drop = TRUE) +
    scale_color_manual(values = colorvect,
                       drop = TRUE)

  #limit y axis scaling to go only 10x smaller than smallest observation
  new_y_min <- fitsub[, 0.1*min(ConcDose,
                            na.rm = TRUE)]
  if(is.finite(new_y_min)){
    p_all <- p_all +
      scale_y_log10(limits = c(new_y_min, NA))
  }

  #axis labeling
  p_all <- p_all +
    xlab("Time, hr") +
    ylab("Concentration/Dose")

  #facet by medium
  p_all <- p_all +
    facet_grid(rows = vars(Media),
               labeller = label_both)

  #and add title
  p_all <- p_all +
    ggtitle(plot_title)

if(save_plot %in% TRUE){

  #and save
  do.call(ggasve,
          args = c(list(filename = filename,
                        plot = p_all,
                        device = file_format,
                        path = file_path),
                   ggsave_args)
          )
}

  if(return_plot %in% TRUE){
    return(p_all)
  }else{
    return(0)
  }

}
