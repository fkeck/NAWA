

#### ML SETUP ####

rf_mod_status <-
  rand_forest(trees = 1000) %>%
  set_engine("ranger") %>%
  set_args(importance = "impurity") %>%
  set_mode("classification")

rf_wf_status <-
  workflow() %>%
  add_model(rf_mod_status) %>%
  add_formula(status ~ .)

null_mod_status <-
  null_model(mode = "classification") %>%
  set_engine("parsnip")

null_wf_status <-
  workflow() %>%
  add_model(null_mod_status) %>%
  add_formula(status ~ .)

# Prepare CDMs

if(project == "Trend") {
  kicknet_IBCH_LU_CDM <- KN_cdm_IBCH %>%
    filter(project == "Trend") %>%
    group_by(sample) %>%
    mutate(indiv_prop = count/sum(count)) %>%
    ungroup() %>%
    rename(site_code = sample)

  edna_LU_CDM <- samples_meta %>%
    filter(project == "NAWA-TREND",
           type == "sample") %>%
    left_join(cdm_OTU97) %>%
    filter(!is.na(count)) %>%
    rename(DNA_seq = OTU) %>%
    group_by(site_code, DNA_seq) %>%
    summarise(count = sum(count)) %>%
    filter(count != 0) %>%
    ungroup() %>%
    flexitarian::rarefy_long(site_code, DNA_seq, count, 10^4) %>%
    group_by(site_code) %>%
    mutate(seq_prop = count/sum(count)) %>%
    ungroup()
}


common_sites <- intersect(unique(kicknet_IBCH_LU_CDM$site_code),
                          unique(edna_LU_CDM$site_code))


# If only Rhine Catchment
if (rhine_catchment_only) {
  common_sites <- intersect(common_sites, subset_sites_rhein$site_code)
}


#### LANDUSE ####

landuse_LU <- landuse_all %>%
  filter(layer == "LU18_10") %>%
  mutate(description = case_when(description %in% c("Aires de bâtiments",
                                                    "Surfaces de transport",
                                                    "Surfaces d'infrastructure spéciale") ~ "Urban",
                                 description %in% c("Arboriculture, viticulture, horticulture",
                                                    "Cultures fourragères et de plein champ") ~ "Agriculture",
                                 description %in% c("Forêt (exploitation agricole non comprise)",
                                                    "Alpages") ~ "Natural")) %>%
  filter(!is.na(description)) %>%
  mutate(description = paste0("fraction_", description)) %>%
  group_by(layer, site_code, description) %>%
  summarise(n = sum(n)) %>%
  group_by(layer, site_code) %>%
  mutate(lu_prop = n/sum(n)) %>%
  ungroup() %>%
  select(site_code, description, lu_prop) %>%
  pivot_wider(names_from = description, values_from = lu_prop) %>%
  filter(site_code %in% common_sites) %>%
  mutate(across(starts_with("fraction_"), ~ ifelse(is.na(.x), 0, .x)))


# FRACTION NATURAL
if(fraction_type == "Natural") {
  landuse_LU <-
    landuse_LU %>%
    mutate(status_q50 = case_when(fraction_Natural < median(landuse_LU$fraction_Natural) ~ "Impacted",
                                  fraction_Natural > median(landuse_LU$fraction_Natural) ~ "Reference")) %>%
    select(site_code, starts_with("status_"))
}

###

landuse_frac <- names(landuse_LU)[-1]


#### KICKNET IBCH ####

KN_IBCH_LU <- kicknet_IBCH_LU_CDM %>%
  dplyr::select(site_code, taxon, indiv_prop) %>%
  pivot_wider(names_from = taxon, values_from = indiv_prop, values_fill = 0) %>%
  inner_join(landuse_LU, .)

KN_IBCH_LU <-
  tibble(method = "KN_IBCH",
         predicted = landuse_frac,
         data = map(landuse_frac, function(x)
           KN_IBCH_LU %>%
             filter(site_code %in% common_sites) %>%
             dplyr::select(site_code, all_of(c(x, names(KN_IBCH_LU)[!str_starts(names(KN_IBCH_LU), "site_code|status_")]))) %>%
             rename_with(function(x) {"status"}, starts_with("status_")) %>%
             filter(!is.na(status)) %>%
             mutate(status = as.factor(status))
         )) %>%
  mutate(folds = map(data, function(x) x %>%
                       dplyr::select(-site_code) %>%
                       vfold_cv(v = 10, repeats = 5)
  ))

KN_IBCH_LU <-
  KN_IBCH_LU %>%
  mutate(rf_fit = map(folds, fit_resamples,
                      object = rf_wf_status,
                      metrics = metric_set(accuracy),
                      control = control_resamples(save_pred = TRUE,
                                                  extract = function(x) ranger::importance(extract_fit_engine(x)))),
         rf_fit_rand = map(folds, fit_resamples,
                           object = null_wf_status,
                           metrics = metric_set(accuracy),
                           control = control_resamples(save_pred = TRUE,
                                                       extract = function(x) ranger::importance(extract_fit_engine(x)))
         ))


KN_IBCH_LU <- bind_rows(
  dplyr::select(KN_IBCH_LU, -rf_fit_rand) %>%
    mutate(null_mod = "FALSE", .after = predicted),
  dplyr::select(KN_IBCH_LU, -rf_fit) %>%
    mutate(null_mod = "TRUE", .after = predicted) %>%
    rename(rf_fit = rf_fit_rand)) %>%
  mutate(metrics_summary = map(rf_fit, collect_metrics),
         metrics = map(rf_fit, collect_metrics, summarize = FALSE),
         predictions = map(rf_fit, collect_predictions),
         var_imp = map2(rf_fit, null_mod, function(x, y) {
           if(y) return(NA)
           x$.extracts %>%
             bind_rows() %>%
             .$.extracts %>%
             map(enframe) %>%
             bind_rows() %>%
             group_by(name) %>%
             summarise(mean_varimp = mean(value))
         }))

KN_IBCH_LU <- select(KN_IBCH_LU, -rf_fit, -folds)




#### EDNA ####

edna_LU <- edna_LU_CDM %>%
  dplyr::select(site_code, DNA_seq, seq_prop) %>%
  pivot_wider(names_from = DNA_seq, values_from = seq_prop, values_fill = 0) %>%
  inner_join(landuse_LU, .)

edna_LU <-
  tibble(method = "eDNA",
         predicted = landuse_frac,
         data = map(landuse_frac, function(x)
           edna_LU %>%
             filter(site_code %in% common_sites) %>%
             dplyr::select(site_code, all_of(c(x, names(edna_LU)[!str_starts(names(edna_LU), "site_code|status_")]))) %>%
             rename_with(function(x) {"status"}, starts_with("status_")) %>%
             filter(!is.na(status)) %>%
             mutate(status = as.factor(status))
         )) %>%
  mutate(folds = map(data, function(x) x %>%
                       dplyr::select(-site_code) %>%
                       vfold_cv(v = 10, repeats = 5)
  ))

edna_LU <-
  edna_LU %>%
  mutate(rf_fit = map(folds, fit_resamples,
                      object = rf_wf_status,
                      metrics = metric_set(accuracy),
                      control = control_resamples(save_pred = TRUE,
                                                  extract = function(x) ranger::importance(extract_fit_engine(x)))),
         rf_fit_rand = map(folds, fit_resamples,
                           object = null_wf_status,
                           metrics = metric_set(accuracy),
                           control = control_resamples(save_pred = TRUE,
                                                       extract = function(x) ranger::importance(extract_fit_engine(x)))
         ))


edna_LU <- bind_rows(
  dplyr::select(edna_LU, -rf_fit_rand) %>%
    mutate(null_mod = "FALSE", .after = predicted),
  dplyr::select(edna_LU, -rf_fit) %>%
    mutate(null_mod = "TRUE", .after = predicted) %>%
    rename(rf_fit = rf_fit_rand)) %>%
  mutate(metrics_summary = map(rf_fit, collect_metrics),
         metrics = map(rf_fit, collect_metrics, summarize = FALSE),
         predictions = map(rf_fit, collect_predictions),
         var_imp = map2(rf_fit, null_mod, function(x, y) {
           if(y) return(NA)
           x$.extracts %>%
             bind_rows() %>%
             .$.extracts %>%
             map(enframe) %>%
             bind_rows() %>%
             group_by(name) %>%
             summarise(mean_varimp = mean(value))
         }))

edna_LU <- select(edna_LU, -rf_fit, -folds)


# Bind all together and produce plots
all_dat_LU <- bind_rows(KN_IBCH_LU, edna_LU)
res_ML[[fraction_type]][["data_results"]] <- all_dat_LU


# Accuracy
sample_location <- landuse_LU %>%
  left_join(sites_meta) %>%
  sf::st_as_sf(coords = c("x_lv03", "y_lv03"), crs = 21781)

sample_success <- all_dat_LU %>%
  filter(null_mod == FALSE, predicted == status) %>%
  select(method, predicted, data, predictions) %>%
  unnest(predictions) %>%
  mutate(site_code = map2_chr(data, .row, function(x, y) {
    x[y, ]$site_code
  })) %>%
  mutate(success = .pred_class == status) %>%
  group_by(site_code, method) %>%
  summarise(success_rate = sum(success)/length(success))

sample_location <- sample_location %>%
  left_join(sample_success)
