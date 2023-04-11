
# Read track
read_csv("data/Trend_track_reads.csv") %>%
  left_join(samples_meta, by = c("Sample" = "sample")) %>%
  filter(site_code %in% common_sites) %>%
  .$Raw %>% sum()

read_csv("data/Trend_track_reads.csv") %>%
  left_join(samples_meta, by = c("Sample" = "sample")) %>%
  filter(site_code %in% common_sites) %>%
  group_by(site_code) %>%
  summarise(Raw = sum(Raw)) %>%
  .$Raw %>% summary()

# Kicknet

res_KN <- KN_IBCH_LU$data[[1]] %>%
  select(-status) %>%
  flexitarian::tbl_to_df(row.names = site_code) %>%
  flexitarian::remove_zerosum_cols()

# Gamma diversity
ncol(res_KN)

# Alpha diversity
res_KN %>% vegan::specnumber() %>% summary()

res_KN %>%
  vegan::specnumber() %>%
  enframe() %>%
  left_join(KN_IBCH_LU$data[[1]], by = c("name" = "site_code")) %>%
  group_by(status) %>%
  summarise(mean(value))





# EDNA

res_edna <- edna_LU$data[[1]] %>%
  select(-status) %>%
  flexitarian::tbl_to_df(row.names = site_code) %>%
  flexitarian::remove_zerosum_cols()

# Gamma diversity
ncol(res_edna)

cdm_OTU97 %>%
  ungroup() %>%
  distinct(OTU, clustered_seq) %>%
  filter(OTU %in% colnames(res_edna)) %>%
  mutate(n_ASV = map_int(clustered_seq, length)) %>%
  summarise(sum(n_ASV))


# Alpha diversity
res_edna %>% vegan::specnumber() %>% summary()

res_edna %>%
  vegan::specnumber() %>%
  enframe() %>%
  left_join(KN_IBCH_LU$data[[1]], by = c("name" = "site_code")) %>%
  group_by(status) %>%
  summarise(mean(value))

# Taxonomic affiliation

cdm_OTU97 %>%
  ungroup() %>%
  distinct(OTU, clustered_seq) %>%
  .$clustered_seq %>%
  Reduce(c, .) %>%
  enframe(value = "DNA_seq") %>%
  left_join(taxa_id) %>%
  apply(2, function(x) sum(!is.na(x)))

cdm_OTU97_expanded %>%
  filter(OTU %in% colnames(res_edna)) %>%
  apply(2, function(x) sum(!is.na(x)))

cdm_OTU97_expanded %>%
  filter(OTU %in% colnames(res_edna)) %>%
  apply(2, function(x) sum(!is.na(x))) %>%
  magrittr::divide_by(max(.))



## T-test accuracy comparisons
accuracies_wide <- all_dat_LU %>%
  select(method, predicted, null_mod, metrics) %>%
  unnest(metrics) %>%
  pivot_wider(id_cols = id:id2, names_from = method:null_mod, values_from = ".estimate")

accuracies_wide %>%
  select(eDNA_status_q50_FALSE,
         KN_IBCH_status_q50_FALSE,
         KN_IBCH_status_q50_TRUE) %>%
  summary()

t.test(accuracies_wide$eDNA_status_q50_FALSE, accuracies_wide$KN_IBCH_status_q50_FALSE, paired = TRUE)
t.test(accuracies_wide$eDNA_status_q50_FALSE, accuracies_wide$KN_IBCH_status_q50_TRUE, paired = TRUE)
t.test(accuracies_wide$KN_IBCH_status_q50_FALSE, accuracies_wide$KN_IBCH_status_q50_TRUE, paired = TRUE)
