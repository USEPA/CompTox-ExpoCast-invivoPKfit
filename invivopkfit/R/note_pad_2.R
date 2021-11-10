# # one_comp <- read.csv("data/PK.fit.table.1comp.alpha.fix.csv")
# # two_comp <- read.csv("data/PK.fit.table.2comp.new.alpha.fix.csv")
#
# one_comp <- PK.fit.table.1comp.new[[1]]
# two_comp <- PK.fit.table.2comp.new[[1]]
# #
# one_comp <- one_comp %>%
#   filter(AIC == Inf | AIC == 200008)
# one_comp <- one_comp %>% distinct(CAS, Reference, Compound, LogLikelihood, AIC)
#
# one_comp_big_sigma <- PK.fit.table.1comp.new[[1]]
# one_comp_big_sigma <- one_comp_big_sigma %>%
#   filter(AIC == Inf | AIC == 200008)
# one_comp_big_sigma <- one_comp_big_sigma %>% distinct(CAS, Reference, Compound, LogLikelihood, AIC)
#
#
# two_comp <- two_comp %>%
#   filter(AIC == Inf | AIC == 200012)
# two_comp <- two_comp %>% distinct(CAS, Reference, Compound, LogLikelihood, AIC)
#
#
#
#
# one_comp_old <- read.csv("data/PK.fit.table.1comp.alpha.fix.csv")
# one_comp_old <- one_comp_old %>%
#   filter(AIC == Inf
#          |
#            LogLikelihood == -99999
#   )
#
# two_comp_old <- read.csv("data/PK.fit.table.2comp.new.alpha.fix.csv")
# two_comp_old <- two_comp_old %>%
#   filter(AIC == Inf
#          |
#            LogLikelihood == -99999
#   )
#
# in_both <- intersect(one_comp$Compound, two_comp$Compound)
#
# # all <- union(one_comp$Compound, two_comp$Compound)
# #
# #
# # benzo_p <- ggplot(data = PK.fit.table.1comp.new[[2]]) +
# #   ggtitle("benzo(a)pyrene, 50-32-8") +
# #   geom_point(data = PK.fit.table.1comp.new[[2]],
# #              aes(x = Time, y = Value))
# # benzo_p <- benzo_p + scale_y_continuous(trans='log10')
# #
# # ondan_p <- ggplot(data = PK.fit.table.1comp.new[[2]]) +
# #   ggtitle("ondansetron, 99614-02-5") +
# #   geom_point(data = PK.fit.table.1comp.new[[2]],
# #              aes(x = Time, y = Value))
# # ondan_p <- ondan_p + scale_y_continuous(trans='log10')
# #
# # bensul_p <- ggplot(data = PK.fit.table.1comp.new[[2]]) +
# #   ggtitle("bensulide, 741-58-2") +
# #   geom_point(data = PK.fit.table.1comp.new[[2]],
# #              aes(x = Time, y = Value))
# # bensul_p <- bensul_p + scale_y_continuous(trans='log10')
# #
# # phenol_p <- ggplot(data = PK.fit.table.1comp.new[[2]]) +
# #   ggtitle("phenolphthalein, 77-09-8") +
# #   geom_point(data = PK.fit.table.1comp.new[[2]],
# #              aes(x = Time, y = Value))
# # phenol_p <- phenol_p + scale_y_continuous(trans='log10')
# #
# # pentachloro_p <- ggplot(data = PK.fit.table.1comp.new[[2]]) +
# #   ggtitle("3,3',4,4',5-pentachlorobiphenyl, 57465-28-8") +
# #   geom_point(data = PK.fit.table.1comp.new[[2]],
# #              aes(x = Time, y = Value))
# # pentachloro_p <- pentachloro_p + scale_y_continuous(trans='log10')
# #
# # pdf("data/not_fit_compounds.pdf")
# # pdf.options(width = 9, height = 7)
# # for (i in 1:length(all_plots)){
# #   print(all_plots[[i]])
# # }
# dev.off()

make_df <- function(data) {


  ### import each sheet from workbook into one list object
  sheet_names <- readxl::excel_sheets(data)

  sheet_list <- lapply(sheet_names,
                       function(i)
                         readxl::read_excel(data, sheet = i)
  )

  names(sheet_list) <- sheet_names

  ### add 'id' columns for 'Documents' sheet and 'Studies' sheet to join by,
  ### as new template format does away with these columns
  sheet_list$Studies$fk_extraction_document_id <- 1

  ### rename 'id' and corresponding 'fk_blank_id' columns to matching name
  sheet_list$Documents <- dplyr::rename(sheet_list$Documents, document_id = id)

  sheet_list$Studies <- dplyr::rename(sheet_list$Studies, study_id = id)

  sheet_list$Studies <- dplyr::rename(sheet_list$Studies, document_id = fk_extraction_document_id)

  sheet_list$Subjects <- dplyr::rename(sheet_list$Subjects, subject_id = id)

  sheet_list$Series <- dplyr::rename(sheet_list$Series, series_id = id)

  sheet_list$Series <- dplyr::rename(sheet_list$Series, study_id = fk_study_id)

  sheet_list$Series <- dplyr::rename(sheet_list$Series, subject_id = fk_subject_id)

  sheet_list$Conc_Time_Values <- dplyr::rename(sheet_list$Conc_Time_Values, series_id = fk_series_id)

  ### combine the five sheets into one data frame using renamed 'id' columns
  df <- sheet_list$Conc_Time_Values %>%
    dplyr::left_join(sheet_list$Series, by = "series_id") %>%
    dplyr::left_join(sheet_list$Studies, by = "study_id") %>%
    dplyr::left_join(sheet_list$Documents, by = "document_id") %>%
    dplyr::left_join(sheet_list$Subjects, by ="subject_id")

  return(df)
}

pmid.8132164 <- make_df("AK submissions/CvT_data_8132164_articles_AK.xlsx")
pmid.32305862 <- make_df("AK submissions/CvT_data_32305862_articles_AK.xlsx")
pmid.32892662 <- make_df("AK submissions/CvT_data_32892662_articles_AK.xlsx")

correct_df <- function(data) {
  data$fk_extraction_document_id <- 1
  data$dsstox_substance_id <- 1

  data <- subset(data, select = c(fk_extraction_document_id,
                                                  dsstox_substance_id,
                                                  series_id,
                                                  time,
                                                  conc,
                                                  analyte_name,
                                                  analyte_casrn,
                                                  time_units,
                                                  conc_units,
                                                  study_id,
                                                  fk_reference_document_id,
                                                  dose_level,
                                                  dose_level_units,
                                                  administration_route,
                                                  conc_medium,
                                                  weight,
                                                  weight_units,
                                                  species))

  data <- data[, c("analyte_casrn",
                   "analyte_name",
                   "dsstox_substance_id",
                   "series_id",
                   "study_id",
                   "fk_reference_document_id",
                   "fk_extraction_document_id",
                   "species",
                   "weight",
                   "weight_units",
                   "dose_level",
                   "dose_level_units",
                   "administration_route",
                   "conc_medium",
                   "conc",
                   "conc_units",
                   "time",
                   "time_units")]

  data <- data %>%
    rename(dsstox_casrn = analyte_casrn,
           analyte_name_original = analyte_name,
           fk_series_id = series_id,
           fk_study_id = study_id,
           weight_kg = weight,
           dose_level_normalized = dose_level,
           administration_route_normalized = administration_route,
           conc_medium_normalized = conc_medium,
           time_hr = time)

  return(data)
}

pmid.8132164 <- correct_df(pmid.8132164)
pmid.8132164$dsstox_casrn <- "51-03-6"
pmid.8132164$dsstox_substance_id <- "DTXSID1021166"
pmid.8132164$conc_medium_normalized <- "blood"
pmid.8132164$fk_reference_document_id[is.na(pmid.8132164$fk_reference_document_id)] <- pmid.8132164$fk_extraction_document_id[is.na(pmid.8132164$fk_reference_document_id)]


pmid.32305862 <- correct_df(pmid.32305862)
pmid.32305862$dsstox_casrn <- "1071-83-6"
pmid.32305862$dsstox_substance_id <- "DTXSID1024122"
pmid.32305862$conc_medium_normalized <- "blood"
pmid.32305862$fk_reference_document_id[is.na(pmid.32305862$fk_reference_document_id)] <- pmid.32305862$fk_extraction_document_id[is.na(pmid.32305862$fk_reference_document_id)]


pmid.32892662 <- correct_df(pmid.32892662)
pmid.32892662$dsstox_casrn <- "1071-83-6"
pmid.32892662$dsstox_substance_id <- "DTXSID1024122"
pmid.32892662$fk_reference_document_id[is.na(pmid.32892662$fk_reference_document_id)] <- pmid.32892662$fk_extraction_document_id[is.na(pmid.32892662$fk_reference_document_id)]



# files <- list.files(path="AK submissions/1_qa_format_complete", pattern="*.xlsx", full.names=TRUE, recursive=FALSE)
#
# big_df <- lapply(files, function(x) {
#   read_excel(x, sheet = "Series")
# })
#
# big_df <- do.call(rbind, big_df)



