
# Figure 1
main_edna_taxa <- c("Ephemeroptera", "Plecoptera", "Coleoptera", "Diptera", "Unclassified", "Calanoida", "Vannellidae undef. order")

treemap_edna <- res_edna %>%
  colSums() %>%
  enframe(name = "OTU", value = "n") %>%
  mutate(n = round(n * 1/min(n))) %>%
  mutate(nk = map2(OTU, n, function(x, y) rep(x, each = y))) %>%
  unnest(nk) %>%
  left_join(cdm_OTU97_expanded) %>%
  replace_na(list(order = "Unclassified")) %>%
  mutate(order = ifelse(order %in% main_edna_taxa, order, "Other")) %>%
  refdb::refdb_set_fields(taxonomy = c("superkingdom" = "superkingdom",
                                       "kingdom" = "kingdom",
                                       "phylum" = "phylum",
                                       "class" = "class",
                                       "order" = "order",
                                       "family" = "family",
                                       "genus" ="genus",
                                       "species" = "OTU")) %>%
  refdb::refdb_plot_tax_treemap(cols = c("order", "OTU"), freq_labels = c(0.01, 1))



main_KN_taxa <- c("Ephemeroptera", "Gastropoda", "Oligochaeta", "Trichoptera", "Plecoptera", "Coleoptera", "Amphipoda", "Diptera")

treemap_KN <- res_KN %>%
  colSums() %>%
  enframe(name = "taxon", value = "n") %>%
  mutate(n = round(n * 1/min(n))) %>%
  mutate(nk = map2(taxon, n, function(x, y) rep(x, each = y))) %>%
  unnest(nk) %>%
  left_join(KN_cdm_IBCH %>% distinct(taxon, .keep_all = TRUE)) %>%
  mutate(group_2 = ifelse(group_2 %in% main_KN_taxa, group_2, "Other")) %>%
  refdb::refdb_set_fields(taxonomy = c("class" = "group_1",
                                       "order" = "group_2",
                                       "family" = "taxon")) %>%
  refdb::refdb_plot_tax_treemap(cols = c("group_2", "taxon"))

library(patchwork)
treemap_KN + treemap_edna


# Figure 2 (4 x 4)
all_dat_LU %>%
  filter(predicted == status) %>%
  mutate(model = c("Morpho-taxonomy", "Null", "eDNA", "Exclude")) %>%
  filter(model != "Exclude") %>%
  unnest(metrics_summary) %>%
  ggplot(aes(x = model, y = mean, ymin = mean - std_err, ymax = mean + std_err)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(position = position_dodge(0.9), width = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Model") +
  ylab("Average accuracy") +
  theme_bw()


# Figure 3 (5*7)
sample_location_bi <-
  sample_location %>%
  select(site_code, status_q50, method, success_rate) %>%
  pivot_wider(names_from = method, values_from = success_rate) %>%
  biscale::bi_class(KN_IBCH, eDNA, style = "equal", dim = 4) %>%
  sf::st_set_geometry("geometry")

rhein_catchment <- sf::read_sf("/home/ecoadmin/Documents/postdoc/Swiss_catchments/EZG_Gewaesser.gpkg", layer = "Teileinzugsgebiet") %>%
  sf::st_transform(21781) %>%
  filter(FLUSSGB %in% c("Rhein", "Aare", "Reuss", "Limmat")) %>%
  sf::st_union(is_coverage = TRUE) %>%
  sf::st_intersection(sf::read_sf("data/Switzerland_shapefile/swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET.shp")) %>%
  sf::st_union(is_coverage = TRUE)

CH_rivers <-
  c("data/Switzerland_shapefile/K4_flusyyyymmdd/k4flusyyyymmdd11_ch2007.shp",
    "data/Switzerland_shapefile/K4_flusyyyymmdd/k4flusyyyymmdd22_ch2007.shp",
    "data/Switzerland_shapefile/K4_flusyyyymmdd/k4flusyyyymmdd33_ch2007.shp",
    "data/Switzerland_shapefile/K4_flusyyyymmdd/k4flusyyyymmdd44_ch2007.shp",
    "data/Switzerland_shapefile/K4_flusyyyymmdd/k4flusyyyymmdd55_ch2007.shp") %>%
  map(sf::read_sf) %>%
  bind_rows() %>%
  sf::st_transform(21781) %>%
  sf::st_intersection(rhein_catchment)

CH_lakes <-
  c("data/Switzerland_shapefile/K4_seenyyyymmdd/k4seenyyyymmdd11_ch2007Poly.shp",
    "data/Switzerland_shapefile/K4_seenyyyymmdd/k4seenyyyymmdd22_ch2007Poly.shp") %>%
  map(sf::read_sf) %>%
  bind_rows() %>%
  sf::st_transform(21781) %>%
  sf::st_intersection(rhein_catchment)


map <- sf::read_sf("data/Switzerland_shapefile/swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET.shp") %>%
  ggplot() +
  geom_sf(color = NA, fill = "grey") +
  geom_sf(fill = NA, size = 0.2, data = rhein_catchment) +
  geom_sf(data = CH_rivers, color = "lightblue", size = 0.3) +
  geom_sf(data = CH_lakes, fill = "lightblue", size = NA) +
  geom_sf(aes(color = bi_class, shape = status_q50), data = sample_location_bi, size = 2.5) +
  biscale::bi_scale_color(pal = "GrPink2", dim = 4) +
  theme_bw() +
  ggspatial::annotation_scale(location = "br", width_hint = 0.15) +
  ggspatial::annotation_north_arrow(location = "tl",   width = unit(1, "cm")) +
  labs(shape = "Status") +
  guides(color = "none")

legend <- biscale::bi_legend(pal = "GrPink2",
                    dim = 4,
                    xlab = "Morphotaxonomy",
                    ylab = "eDNA",
                    size = 8, arrows = FALSE,
                    breaks = list(bi_x = c("0%", "50%", "", "100%"), bi_y = c("0%", "50%", "", "100%")))
cowplot::ggdraw() +
  cowplot::draw_plot(map, 0, 0, 1, 1) +
  cowplot::draw_plot(legend, 0.79, .12, 0.25, 0.25)

# Figure 4
all_dat_LU %>%
  filter(null_mod == FALSE, predicted == status) %>%
  select(method, var_imp) %>%
  mutate(method = c("Morpho-taxonomy", "eDNA")) %>%
  unnest(var_imp) %>%
  group_by(method) %>%
  slice_max(order_by = mean_varimp, n = 10) %>%
  mutate(name = str_remove_all(name, "`")) %>%
  left_join(cdm_OTU97_expanded, by = c("name" = "OTU")) %>%
  left_join(select(KN_cdm_IBCH, group_1, group_2, taxon) %>% distinct(), by = c("name" = "taxon")) %>%
  mutate(lowest_taxa = ifelse(str_starts(lowest_taxa, "[0-9]"), "Unclassified", lowest_taxa)) %>%
  mutate(lowest_taxa = ifelse(lowest_level > 5, paste0(order, ": ", lowest_taxa), lowest_taxa)) %>%
  mutate(lowest_taxa = ifelse(!is.na(lowest_level), paste0(name, " (", lowest_taxa, ")"), lowest_taxa)) %>%
  mutate(lowest_taxa = ifelse(is.na(lowest_level), paste0(group_2, ": ", name), lowest_taxa)) %>%
  ggplot() +
  geom_col(aes(fct_reorder(lowest_taxa, mean_varimp), mean_varimp)) +
  facet_wrap(vars(method), scales = "free", ncol = 2) +
  coord_flip() +
  xlab("") +
  ylab("Variable Importance") +
  theme(axis.title = element_blank()) +
  theme_bw()



# Figure SI-1

landuse_all %>%
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
  mutate(across(starts_with("fraction_"), ~ ifelse(is.na(.x), 0, .x))) %>%
  ggplot() +
  geom_histogram(aes(fraction_Natural), bins = 20) +
  geom_vline(aes(xintercept = median(.data$fraction_Natural))) +
  scale_x_continuous(labels = scales::percent) +
  theme_bw() +
  xlab("Proportion of natural land use") +
  ylab("Number of sites")
