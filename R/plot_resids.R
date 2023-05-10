#' Plot residuals from winning model
#'
#' Make diagnostic plots of residuals from the winning TK model for a dataset
#'
#' For a given dataset (DTXSID and Species), this function calculates the
#' predicted concentrations from the winning (best-fit by lowest AIC) TK model.
#' Then it computes the residuals (observed - predicted), both dose-normalized
#' and non-dose-normalized. Then it creates a grid of plots with two rows and
#' three columns. The first row shows non-dose-normalized residuals; the second
#' row shows dose-normalized residuals. The first column plots residuals vs.
#' observed concentration. The second column plots residuals vs. dose. The third
#' column plots residuals vs. time.
#'
#' These plots can be used to diagnose heteroscedasticity in residual errors
#' from the best-fit model.
#'
#' @param cvt_pre Original concentration vs. time data, pre-processed using
#'   `preprocess_data()`.
#' @param fit_flat Fitting output for the flat model: the output of `fit_all()`
#'   with `model = "flat"`
#' @param fit_1comp Fitting output for the 1-compartment model: the output of
#'   `fit_all()` with `model = "1compartment"`
#' @param fit_2comp Fitting output for the 2-compartment model: the output of
#'   `fit_all()` with `model = "2compartment"`
#' @return A `data.table` of goodness-of-fit measures for each dataset, each
#'   model, and each analysis (joint, separate, and pooled)
#' @author Caroline Ring
#' @export
#' @import data.table patchwork
plot_resids <- function(cvt_pre,
                         fit_flat,
                         fit_1comp,
                         fit_2comp,
                        return_plots = FALSE,
                        save_indiv_plots = TRUE,
                        save_combined_plots = TRUE,
                        file_path = "../inst/ext/plots/resid_plots",
                        file_suffix = NULL,
                        file_format = "pdf",
                        ggsave_args = list(scale = 1,
                                           width = 13.33,
                                           height = 7.5,
                                           units = "in",
                                           dpi = 300,
                                           limitsize = TRUE,
                                           bg = NULL)){

 pred_DT <- get_all_predictions(cvt_pre =cvt_pre,
                               fit_flat = fit_flat,
                               fit_1comp = fit_1comp,
                               fit_2comp = fit_2comp)

 #Get winning model
 pk_fit <- merge_fits(fit_flat = fit_flat,
                      fit_1comp = fit_1comp,
                      fit_2comp = fit_2comp)
 pk_win <- unique(pk_fit[Analysis_Type  %in% "Joint",
                         .(DTXSID, Species, winmodel)])
 setnames(pk_win, "winmodel", "model")


 #Keep only the predictions that match the winning models
 pred_DT <- pk_win[pred_DT, on = names(pk_win), nomatch = 0]


  #What do residuals look like when plotted vs. concentration? Vs.
  #dose-normalized concentration? Vs dose? Vs time? Is there heteroscedasticity?

  pred_DT[, resid := Conc - pred_conc]
  pred_DT[, resid_Dose := Conc_Dose - pred_conc/Dose]

  pred_DT <- pred_DT[Detect %in% "Detect"]

  tmp <- unique(pred_DT[Analysis_Type %in% "Joint",
                        .(DTXSID, Species, Analysis_Type, Studies.Analyzed, model)])

  setnames(tmp,
           names(tmp),
           paste(names(tmp),
                 "in",
                 sep = "_"))

plot_list <- mapply(function(DTXSID_in,
                             Species_in,
                             Analysis_Type_in,
                             Studies.Analyzed_in,
                             model_in){
  dat <- pred_DT[DTXSID %in% DTXSID_in &
                   Species %in% Species_in &
                   Analysis_Type %in% Analysis_Type_in &
                   Studies.Analyzed %in% Studies.Analyzed_in &
                   model %in% model_in]
  # dat[, Detect:= factor(Detect,
  #                       levels = c("Detect", "Non-Detect"))]
  dat[, Route := factor(Route,
                        levels = c("iv", "po"))]
  dat[, Media := factor(Media, levels = c("blood", "plasma"))]
  dat[, Route_Media := interaction(Route, Media, drop = FALSE)]
  x1vars <- c("Conc", "Dose", "Time")
  x2vars <- c("Conc_Dose", "Dose", "Time")
  row1 <- lapply(x1vars,
                 function(this_xvar){
                   xlabel <- ifelse(this_xvar %in% "Conc",
                          "Conc, mg/L",
                          ifelse(this_xvar %in% "Dose",
                                 "Dose, mg/kg",
                                 "Time, hours"))
                   ggplot(data = dat,
                          aes_string(x = this_xvar,
                              y = "resid",
                              color = "Route_Media",
                              shape = "Reference"
                          )) + geom_blank() +
                     ylab("(Obs - Pred), mg/L") +
                     xlab(xlabel)
                 }
  )

  row2 <- lapply(x2vars,
                        function(this_xvar){
                          xlabel <- ifelse(this_xvar %in% "Conc_Dose",
                                           "Conc/Dose, L/kg",
                                           ifelse(this_xvar %in% "Dose",
                                                  "Dose, mg/kg",
                                                  "Time, hours"))
                          ggplot(data = dat,
                                 aes_string(x = this_xvar,
                                            y = "resid_Dose",
                                            color = "Route_Media",
                                            shape = "Reference"
                                 )) + geom_blank() +
                            ylab("(Obs - Pred)/Dose, kg/L") +
                            xlab(xlabel)
                        }
  )

  all_plots <- c(row1, row2)


all_plots <- lapply(all_plots,
                    function(this_plot){

                    this_plot +
                        geom_point(fill = Route_Media,
                                   shape = 21,
                                   size = 2,
                                   stroke = 1) +
                     # geom_point(fill = "white",
                     #            #color = "gray50",
                     #            shape = 21,
                     #            size = 2,
                     #            stroke = 1) +
                     # geom_point(aes(alpha = Detect,
                     #                fill = Route_Media),
                     #            shape = 21,
                     #            #color = "gray50",
                     #            #fill = "gray50",
                     #            size = 2,
                     #            stroke = 1) +
                        scale_color_brewer(palette = "Paired") +
                        scale_fill_brewer(palette = "Paired") +
                     # scale_alpha_manual(values = c("Detect" = 1, "Non-Detect" = 0),
                     #                    drop = FALSE,
                     #                    name = NULL) +
                     # guides(alpha = guide_legend(override.aes = list(shape = 21,
                     #                                                 color = "black",
                     #                                                 stroke = 1,
                     #                                                 fill = c("black", NA),
                     #                                                 alpha = 1)
                     # )
                     # ) +
                        scale_shape_manual(values = 21:25)  #use only the 5 shapes where a fill can be added
                        #this limits us to visualizing only 5 References
                        #but admittedly it's hard to distinguish more than 5 shapes anyway

                 })

#generate plot title
plot_title <- paste0(DTXSID_in,
                     " (", dat[, unique(Compound)], ")\n",
                     "Species = ", Species_in, ", ",
                     "Doses = ", paste(signif(
                       sort(unique(dat$Dose)),
                       3),
                       collapse = ", "), " mg/kg\n",
                     "Analysis Type = ", paste(Analysis_Type_in, collapse= ", "))

plot_subtitle <- paste0("Winning Model = ",
                        model_in)
#use patchwork to produce grid of plots with overall title & subtitle
p_grid <- wrap_plots(all_plots,
                     ncol = 3,
                     guides = "collect") +
    plot_annotation(title = plot_title,
                    subtitle = plot_subtitle)

if(save_indiv_plots %in% TRUE){
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

  do.call(ggsave,
          args = c(list(filename = filename,
                        plot = p_grid,
                        device = file_format,
                        path = file_path),
                   ggsave_args)
  )
}

if(save_combined_plots %in% TRUE |
   return_plots %in% TRUE){
  return(p_grid)
}else{
  return(0)
}
},
DTXSID = tmp$DTXSID_in,
Species_in = tmp$Species_in,
Analysis_Type_in = tmp$Analysis_Type_in,
Studies.Analyzed_in = tmp$Studies.Analyzed_in,
model_in = tmp$model_in,
SIMPLIFY = FALSE)

if(save_combined_plots %in% TRUE){
  #set up filename for saving
  filename <- "all_resid_plots"

  if(length(file_suffix) > 0){
    filename <- paste0(filename,
                       "_",
                       file_suffix)
  }

  filename <- paste0(filename, ".", file_format)

  pdf(paste(file_path, filename, sep = "/"),
      height = ggsave_args$height,
      width = ggsave_args$width)
  print(plot_list)
  dev.off()
}

if(return_plots %in% TRUE){
  return(plot_list)
}else{
  return(0)
}

}
