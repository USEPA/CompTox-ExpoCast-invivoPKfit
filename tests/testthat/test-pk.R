test_that("creating a pk object works for CvT data for DTXSID3061635",
          {
            expect_no_error(my_pk <- pk(
              data = subset(cvt,
                            chemicals_analyzed.dsstox_substance_id %in% "DTXSID3061635"
              )
            )
            )
          }
)

test_that("pk object has the expected default mapping",
          {
            my_pk <- pk(
              data = subset(cvt,
                            chemicals_analyzed.dsstox_substance_id %in% "DTXSID3061635"
              )
            )
            expect_equivalent(my_pk$mapping,
                              ggplot2::aes(Chemical = chemicals_analyzed.dsstox_substance_id,
                                           DTXSID = chemicals_analyzed.dsstox_substance_id,
                                           Chemical_Name = chemicals_analyzed.preferred_name,
                                           CASRN = chemicals_analyzed.dsstox_casrn,
                                           Species = subjects.species,
                                           Reference = as.character(
                                             ifelse(
                                               is.na(
                                                 documents_reference.id
                                               ),
                                               documents_extraction.id,
                                               documents_reference.id
                                             )
                                           ),
                                           Media = series.conc_medium_normalized,
                                           Route = studies.administration_route_normalized,
                                           Dose = studies.dose_level_normalized,
                                           Dose.Units = "mg/kg",
                                           Subject = subjects.id,
                                           N_Subjects =  series.n_subjects_in_series,
                                           Weight = subjects.weight_kg,
                                           Weight.Units = "kg",
                                           Time = conc_time_values.time_hr,
                                           Time.Units = "hours",
                                           Value = conc_time_values.conc,
                                           Value.Units = "mg/L",
                                           LOQ = series.loq_normalized,
                                           Value_SD  = conc_time_values.conc_sd_normalized
                              ))
          })

test_that("creating a pk object with a mapping missing all required variables throws the expected warning",
          {
            expect_warning(
              pk(
            data = subset(cvt,
                          chemicals_analyzed.dsstox_substance_id %in% "DTXSID3061635"
            ),
            mapping = NULL
              )
            )
          }
          )


