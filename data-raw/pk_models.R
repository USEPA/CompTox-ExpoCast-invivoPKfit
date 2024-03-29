devtools::load_all()

#' 1-compartment model object.
model_1comp <- pk_model(name = "model_1comp",
                    params = c("kelim", "Vdist", "Fgutabs", "kgutabs", "Fgutabs_Vdist", "Rblood2plasma"),
                    conc_fun = "cp_1comp",
                    auc_fun = "auc_1comp",
                    params_fun = "get_params_1comp",
                    tkstats_fun = "tkstats_1comp",
                    conc_fun_args = NULL,
                    auc_fun_args = NULL,
                    params_fun_args = NULL,
                    tkstats_fun_args = NULL)


model_2comp <- pk_model(name = "model_2comp",
                    params = c("kelim", "k12", "k21", "V1", "Fgutabs", "kgutabs", "Fgutabs_V1", "Rblood2plasma"),
                    conc_fun = "cp_2comp",
                    auc_fun = "auc_2comp",
                    params_fun = "get_params_2comp",
                    tkstats_fun = "tkstats_2comp",
                    conc_fun_args = NULL,
                    auc_fun_args = NULL,
                    params_fun_args = NULL,
                    tkstats_fun_args = NULL)

model_flat <- pk_model(name = "model_flat",
                    params = c("Vdist", "Fgutabs", "Fgutabs_Vdist", "Rblood2plasma"),
                    conc_fun = "cp_flat",
                    auc_fun = "auc_flat",
                    params_fun = "get_params_flat",
                   tkstats_fun = "tkstats_flat",
                    conc_fun_args = NULL,
                    auc_fun_args = NULL,
                    params_fun_args = NULL,
                   tkstats_fun_args = NULL)

usethis::use_data(model_1comp,
                  model_2comp,
                  model_flat,
                  internal = FALSE,
                  overwrite = TRUE)

