# qb_df is the CSV file that contains the tested chemicals for the QSAR bakeoff
qb_df <- paste0(Sys.getenv("TKQSAR_DIR"), "TestChemsDashboardInfo.csv")
qb_dtxsid <- qb_df$DTXSID

table(qb_dtxsid %in% cvt_df$analyzed_chem_dtxsid)

qb_df_tracker <- qb_df %>%
  distinct(DTXSID, PREFERRED_NAME, CASRN, SMILES)

qb_df_tracker <- qb_df_tracker %>%
  mutate(include = TRUE,
         reason = "")


rad_experiments <- cvt_df %>%
  filter(radiolabeled, analyzed_chem_dtxsid %in% qb_df$DTXSID) %>%
  pull(analyzed_chem_dtxsid) %>% unique()


qb_df_tracker <- qb_df_tracker %>%
  mutate(
    include = case_match(
      DTXSID,
      rad_experiments ~ FALSE,
      .default = include
    ),
    reason = case_match(
      DTXSID,
      rad_experiments ~ "Radiolabeling experiment.",
      .default = reason
    )
  )

incongruent_chems <-  cvt_df %>%
  rowwise() %>%
  filter(
    !(fk_dosed_chemical_id == fk_analyzed_chemical_id),
    analyzed_chem_dtxsid %in% qb_df$DTXSID
  ) %>%
  pull(analyzed_chem_dtxsid) %>% unique()

congruent_chems <- cvt_df %>%
  rowwise() %>%
  filter(
    identical(fk_dosed_chemical_id, fk_analyzed_chemical_id),
    analyzed_chem_dtxsid %in% qb_df$DTXSID
  ) %>%
  pull(analyzed_chem_dtxsid) %>% unique()

incongruent_chems <- setdiff(incongruent_chems, congruent_chems)


qb_df_tracker <- qb_df_tracker %>%
  mutate(
    include = case_match(
      DTXSID,
      incongruent_chems ~ FALSE,
      .default = include
    ),
    reason = case_match(
      DTXSID,
      incongruent_chems ~ trimws(paste(reason, "Dosed chemical is not analyzed chemical.")),
      .default = reason
    )
  )


qb_df_tracker <- qb_df_tracker %>%
  mutate(
    include = case_match(
      DTXSID,
      qb_dtxsid[!(qb_dtxsid %in% cvt_df$analyzed_chem_dtxsid)] ~ FALSE,
      .default = include
    ),
    reason = case_match(
      DTXSID,
      qb_dtxsid[!(qb_dtxsid %in% cvt_df$analyzed_chem_dtxsid)] ~ trimws(paste(reason, "Missing from CvTdb.")),
      .default = reason
    )
  )


cvt_df_qb <- cvt_df %>%
  filter(analyzed_chem_dtxsid %in% qb_df_tracker[qb_df_tracker$include,][["DTXSID"]])

na_conc_chems <- cvt_df_qb %>%
  filter((is.na(conc) & is.na(conc_original))) %>%
  pull(analyzed_chem_dtxsid) %>% unique()

good_conc_chems <- cvt_df_qb %>%
  filter(!(is.na(conc) & is.na(conc_original))) %>%
  pull(analyzed_chem_dtxsid) %>% unique()

na_conc_chems <- setdiff(na_conc_chems, good_conc_chems)

negative_conc_chems <- cvt_df_qb %>%
  filter(!((conc < 0 | conc_original < 0) %in% FALSE)) %>%
  pull(analyzed_chem_dtxsid) %>% unique()

positive_conc_chems <- cvt_df_qb %>%
  filter(!((conc >= 0 & conc_original >= 0) %in% FALSE)) %>%
  pull(analyzed_chem_dtxsid) %>% unique()

negative_conc_chems <- setdiff(negative_conc_chems, positive_conc_chems)

qb_df_tracker <- qb_df_tracker %>%
  mutate(
    include = case_match(
      DTXSID,
      na_conc_chems ~ FALSE,
      negative_conc_chems ~ FALSE,
      .default = include
    ),
    reason = case_match(
      DTXSID,
      na_conc_chems ~ trimws(paste(reason,
                                   "NA concentration values, with no original concentration values as backup.")),
      negative_conc_chems ~ trimws(paste(reason,
                                         "Negative concentration values.")),
      .default = reason
    )
  )

cvt_df_qb %>%
  filter(!(is.na(conc) & is.na(conc_original))) %>%
  group_by(fk_study_id) %>%
  filter(n() <= 4) %>%
  ungroup() %>%
  pull(analyzed_chem_dtxsid) %>%
  unique() -> low_observation_chem

cvt_df_qb %>%
  filter(!(is.na(conc) & is.na(conc_original)),
         !(conc_medium_normalized %in% c("plasma", "blood"))) %>%
  pull(analyzed_chem_dtxsid) %>%
  unique() -> plbl_chem

plbl_chem <- setdiff(
  plbl_chem,
  subset(cvt_df_qb,
         subset = (conc_medium_normalized %in% c("plasma", "blood")))[["analyzed_chem_dtxsid"]]
  )

qb_df_tracker <- qb_df_tracker %>%
  mutate(
    include = case_match(
      DTXSID,
      low_observation_chem ~ FALSE,
      plbl_chem ~ FALSE,
      .default = include
    ),
    reason = case_match(
      DTXSID,
      low_observation_chem ~ trimws(paste(reason, "Less than 5 observations per study_id.")),
      plbl_chem ~ trimws(paste(reason, "Chemical has no plasma or blood data.")),
      .default = reason
    )
  )



qb_df_tracker %>% count(reason)

## FINISH cvt object creation... run fits. Look at included chemicals in qb_df_tracker
load("data-raw/temp_pk.rda")

table(qb_df_tracker[qb_df_tracker$include, ][["DTXSID"]] %in% cvt_df$analyzed_chem_dtxsid)
qb_df_tracker[qb_df_tracker$include, ][["DTXSID"]] -> included_dtxsid

winmodels <- get_winning_model(my_pk)
count(winmodels, model)

nofit_chems <- setdiff(
  intersect(
    included_dtxsid,
    my_pk$prefit$fit_check %>%
      filter(fit_decision == "abort") %>%
      pull(Chemical)
  ),
  my_pk$prefit$fit_check %>%
    filter(fit_decision != "abort") %>%
    semi_join(winmodels[1:4]) %>%
    pull(Chemical)
)


flat_chems <- setdiff(
  winmodels %>%
    filter(model == "model_flat", Chemical %in% included_dtxsid) %>%
    pull(Chemical),
  winmodels %>%
    filter(model != "model_flat", Chemical %in% included_dtxsid) %>%
    pull(Chemical)
)

qb_df_tracker <- qb_df_tracker %>%
  mutate(
    include = case_match(
      DTXSID,
      nofit_chems ~ FALSE,
      flat_chems ~ FALSE,
      .default = include
    ),
    reason = case_match(
      DTXSID,
      nofit_chems ~ trimws(paste(reason, "Could not fit chemical.")),
      flat_chems ~ trimws(paste(reason, "Null model was best fit.")),
      .default = reason
    )
  )

inf_auc_chems <- setdiff(
  eval_tkstats(my_pk, finite_only = FALSE) %>%
    filter(is.infinite(AUC_infinity.tkstats), Chemical %in% included_dtxsid,
           model != "model_flat") %>%
    pull(Chemical),
  eval_tkstats(my_pk, finite_only = FALSE) %>%
    filter(is.finite(AUC_infinity.tkstats), Chemical %in% included_dtxsid,
           model != "model_flat") %>%
    pull(Chemical)
)


qb_df_tracker <- qb_df_tracker %>%
  mutate(
    include = case_match(
      DTXSID,
      inf_auc_chems ~ FALSE,
      .default = include
    ),
    reason = case_match(
      DTXSID,
      inf_auc_chems ~ trimws(paste(
        reason,
        "In eval_tkstats, AUC_infinity was non-finite."
      )),
      .default = reason
    )
  )

count(qb_df_tracker, reason)

count(qb_df_tracker, include)

included_dtxsid <- unique(qb_df_tracker[qb_df_tracker$include, ][["DTXSID"]])

my_fit <- get_fit(my_pk)

chems_for_qsar_bakeoff <- my_fit %>%
  filter(Chemical %in% included_dtxsid, convcode == 0,
         model %in% c("model_1comp", "model_2comp"),
         Species %in% c("human", "rat")) %>%
  distinct(Chemical, Species) %>%
  semi_join(winmodels)

# Quick aside, which chemicals have both rat and human?
chems_for_qsar_bakeoff %>%
  distinct(Chemical, Species) %>%
  group_by(Chemical) %>%
  count(Chemical, name = "N_Species") %>%
  filter(N_Species == 2)

cvt_qb <- cvt %>%
  semi_join(
    chems_for_qsar_bakeoff %>% distinct(Chemical, Species),
    by = join_by(analyzed_chem_dtxsid == Chemical, species == Species)
  )

# Some quick sanity checks
apply(cvt_qb, 2, \(x) {length(unique(x))})

# There must be at least 5 concentration values above loq
cvt_qb %>% group_by(analyzed_chem_dtxsid, species) %>%
  summarize(below_loq = sum(invivPK_conc <= invivPK_loq),
            above_loq = sum(invivPK_conc > invivPK_loq)) %>%
  arrange(above_loq)

# All looks good. Saving as select CvTdb data
ggplot(cvt_qb,
       aes(x = time_hr, y = invivPK_conc, color = analyzed_chem_dtxsid)) +
  geom_point() +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none")

readr::write_csv(
  x = cvt_qb,
  file = paste0(Sys.getenv("TKQSAR_DIR"), "CvTdb_selectData_2025May.csv")
)

my_tkstats_fin <- eval_tkstats(my_pk, finite_only = FALSE) %>%
  semi_join(chems_for_qsar_bakeoff)

my_fit_fin <- get_fit(my_pk) %>%
  semi_join(chems_for_qsar_bakeoff) %>%
  filter(!startsWith(param_name, "sigma_")) %>%
  select(Chemical:convcode) %>% distinct() %>%
  pivot_wider(id_cols = Chemical:method,
              names_from = param_name,
              values_from = estimate)

my_tkstats_fin <- my_tkstats_fin %>%
  left_join(my_fit_fin)

my_tkstats_fin <- my_tkstats_fin %>%
  select(!ends_with(".nca")) %>% distinct()

my_preds_fin <- predict(my_pk) %>%
  semi_join(chems_for_qsar_bakeoff)

length(unique(my_tkstats_fin$Chemical))
length(unique(my_preds_fin$Chemical))

writexl::write_xlsx(
  x = list(tkstats = my_tkstats_fin,
           predictions = my_preds_fin),
  path = paste0(Sys.getenv("TKQSAR_DIR"), "evalTKstats_bakeoff_2025May.xlsx")
)

# treemap plot of qb_df_tracker

qb_df_tracker %>%
  mutate(reason = ifelse(reason == "", "OK", reason)) %>%
  count(reason, name = "N_Chemicals") %>%
  ggplot(
    aes(
      area = N_Chemicals,
      fill = forcats::fct_relevel(
        paste0(reason, " (", N_Chemicals, ")"),
        "OK (81)"
        ),
      label = N_Chemicals
    )
  ) +
  treemapify::geom_treemap(color = NA) +
  treemapify::geom_treemap_text(color = "grey5", place = "center", grow = TRUE) +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(nrow = 4)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

readr::write_csv(
  x = qb_df_tracker,
  file = paste0(
    Sys.getenv("TKQSAR_DIR"),
    "TestChemsDashboardInfo_FilteringNotes_2025May.csv"
  )
)
