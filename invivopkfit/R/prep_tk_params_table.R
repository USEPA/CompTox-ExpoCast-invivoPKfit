# library(tidyverse)
#
# ### define data
# tk_params <- read.csv("inst/ext/tk_parameters_202202241334.csv")
# pk_data <- pk_output_1comp_04142022
# raw_data <- raw_1comp_04142022
# og_raw_data <- read.csv("inst/ext/cvtdb_invivopkfit_04142022.csv")
#
# ### get rid of predicted parameters
# # will just get in the way
# raw_data <- dplyr::select(raw_data, -c("kelim",
#                                        "Vdist",
#                                        "Fgutabs",
#                                        "kgutabs"))
#
# ### make reference column same class between both data sets before joining
# raw_data$Reference <- as.character(raw_data$Reference)
#
# ### get rid of the stupid row number col
# raw_data$X <- NULL
# pk_data$X <- NULL
#
# pk_raw <- left_join(pk_data, raw_data)
#
# pk_raw <- pk_raw %>%
#   filter(param.value.type == "Fitted geometric mean")
#
# ### right now I'm taking out
# # joint analysis, until I
# # figure out how to better
# # handle it
# pk_raw <- pk_raw %>%
#   filter(Data.Analyzed != "Joint Analysis")
#
# og_raw_data$fk_reference_document_id <- as.character(og_raw_data$fk_reference_document_id)
#
# pk_raw_og <- left_join(pk_raw,
#                        og_raw_data,
#                        by = c("CAS" = "dsstox_casrn",
#                               "Reference" = "fk_reference_document_id",
#                               "Media" = "conc_medium_normalized",
#                               "Species" = "species"))
# pk_raw_og <- pk_raw_og %>%
#   distinct(series_id, .keep_all = TRUE)
#
# pk_raw_og <- pk_raw_og %>% dplyr::select(series_id,
#                                          CAS,
#                                          Species,
#                                          kelim,
#                                          Vdist,
#                                          Fgutabs,
#                                          sigma_value,
#                                          Reference,
#                                          Data.Analyzed,
#                                          Compound,
#                                          LogLikelihood,
#                                          AIC,
#                                          model,
#                                          CLtot,
#                                          Css,
#                                          halflife,
#                                          tpeak.oral,
#                                          Cpeak.oral.1mgkg)
#
# pk_raw_og$id <- seq.int(nrow(pk_raw_og))
#
# pk_raw_og <- pk_raw_og %>%
#   pivot_longer(cols = c("kelim",
#                         "Vdist",
#                         "Fgutabs",
#                         "sigma_value",
#                         "CLtot",
#                         "Css",
#                         "halflife",
#                         "tpeak.oral",
#                         "Cpeak.oral.1mgkg"),
#                names_to = "parameter_name",
#                values_to = "parameter_value")
#
# pk_raw_og <- pk_raw_og %>%
#   rename(fk_series_id = series_id)
#
# pk_raw_og <- pk_raw_og[, c("id",
#                            "CAS",
#                            "Species",
#                            "Reference",
#                            "fk_series_id",
#                            "parameter_name",
#                            "parameter_value",
#                            "LogLikelihood",
#                            "AIC")]
#
# pk_raw_og <- pk_raw_og %>%
#   group_by(CAS, Species, Reference) %>%
#   mutate(series_group_id = cur_group_id())
#
#
# pk_raw_og <- pk_raw_og[, c("id",
#                            "CAS",
#                            "series_group_id",
#                            "fk_series_id",
#                            "parameter_name",
#                            "parameter_value",
#                            "LogLikelihood",
#                            "AIC")]

