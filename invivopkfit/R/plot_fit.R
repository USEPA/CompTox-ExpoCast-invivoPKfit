plot_fit <- function(DTXSID_in,
                     Species_in,
                     Analysis_Type_in,
                     Reference_in,
                     model_in = c("1compartment", "2compartment", "flat"),
                     fitdata,
                     pk_fit){
  #if more than one refrence in comma separated list, get a vector of references
  reflist <- strsplit(Reference_in,
                      split =  ", ")[[1]]

  #subset the data appropriately
  fitsub <- subset(fitdata, DTXSID %in% DTXSID_in &
                     Species %in% Species_in &
                     Reference %in% reflist)

  #subset the fitted parameter data appropriately
  pksub <- subset(pk_fit, DTXSID %in% DTXSID_in &
                    Species %in% Species_in &
                    Reference %in% Reference_in &
                    Analysis_Type %in% Analysis_Type_in &
                    model %in% model_in)

  #Get concentration or LOQ
  fitsub[, Conc:=pmax(Value, LOQ, na.rm = TRUE)]
  #Get detection flag
  fitsub[, Detect:=!is.na(Value)]

  #calculate concentration normalized to dose
  #this allows easier visualization
  fitsub[, ConcDose:=Conc/Dose]

  fitsub[, Dose_cat:=factor(Dose)]

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

  #color will represent dose
  #filled symbols will represent detects (>LOQ)
  #emtpy symbols will represent nondetects (<LOQ)
  #circles represent po, triangles represent IV

  #So we need to map shape to ND
  shapevect <- c(
                 "FALSE" = 1, #open circle
                 "TRUE" = 19 #filled circle
  )


  #plot the conc-time-dose data
  p <- ggplot(data = fitsub,
              aes(x = Time,
                  y = ConcDose,
                  color = Route,
                  shape = Detect
              )) +
    geom_point()

  #get a table of unique combinations of model & route
  modelroute <- expand.grid(model = model_in,
                            route = unique(fitsub$Route),
                            stringsAsFactors = FALSE)

  #For each model & route, add a stat_function() layer to solve the model
  layer_list <- mapply(function(mod, route){
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
                       breaks = c("TRUE.po",
                                  "FALSE.po",
                                  "TRUE.iv",
                                  "FALSE.iv"),
                       labels = c("Oral, nondetect",
                                  "Oral, detect",
                                  "IV, nondetect",
                                  "IV, detect"),
                       drop = TRUE) +
    scale_linetype_manual(name = "Model",
                          values = c("1compartment" = 1,
                                     "2compartment" = 2,
                                     "flat" = 3),
                          drop = TRUE)

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
    ggtitle(paste0(DTXSID_in, "\n",
                   "Species = ", Species_in, "\n",
                   "Reference = ", Reference_in, "\n",
                   "Analysis Type = ", Analysis_Type_in))


  #and save
  ggsave(paste0("inst/ext/plots/",
                DTXSID_in,
                "_",
                Species_in,
                "_",
                Reference_in,
                "_",
                Analysis_Type_in,
                ".pdf"),
         p_all,
         height = 8.5,
         width = 11)

  return(0)

}
