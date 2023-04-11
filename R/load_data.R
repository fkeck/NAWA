library(tidyverse)
library(bioseq)
library(flexitarian)
library(vegan)
library(magrittr)
library(tidymodels)


sites_meta <- read_csv("data/sites_meta.csv")
samples_meta <- read_csv("data/samples_meta.csv") %>%
  filter(project == "NAWA-TREND")
landuse <- read_csv("data/landuse_3cat.csv")
landuse_all <- read_csv("data/landuse.csv")


# Read KN IBCH
cdm_1 <- read_csv("data/kicknet_Trend_IBCH.csv") %>%
  pivot_longer(-1:-3, names_to = "sample", values_to = "count") %>%
  mutate(project = "Trend", .before = group_1)

KN_cdm_IBCH <- cdm_1

#rm(cdm_1, cdm_2, cdm_3)

# Read DNA
taxa_id <- read_csv("data/merged_taxo_IdTaxa.csv")
taxa_conf <- read_csv("data/merged_confidence_IdTaxa.csv")
taxa_id[, 3:10][taxa_conf[, 3:10] < 60] <- NA

taxa_id <- taxa_id %>%
  select(-roottree) %>%
  mutate(DNA_seq = as_dna(DNA_SEQ), DNA_SEQ = NULL, .before = 1)


cdm <- read_csv("data/merged_seqtab_nochim.csv") %>%
  select(DNA_SEQ, any_of(samples_meta$sample)) %>%
  rename(DNA_seq = DNA_SEQ) %>%
  mutate(DNA_seq = as_dna(DNA_seq)) %>%
  pivot_longer(!DNA_seq, names_to = "sample", values_to = "count") %>%
  mutate(count = as.integer(count))

rare_ASV <- cdm %>%
  group_by(DNA_seq) %>%
  summarise(count = sum(count)) %>%
  filter(count <= 10L)

bad_length_ASV <- cdm %>%
  distinct(DNA_seq) %>%
  mutate(read_len = nchar(DNA_seq)) %>%
  filter(read_len < 137 | read_len > 147)

n_control_ASV <- samples_meta %>%
  filter(str_detect(type, "n_control")) %>%
  left_join(cdm) %>%
  group_by(DNA_seq) %>%
  summarise(count = sum(count)) %>%
  left_join(group_by(cdm, DNA_seq) %>% count(wt = count, name = "n")) %>%
  mutate(ratio = count / n) %>%
  filter(ratio > 0.01)

cdm <- cdm %>%
  filter(!DNA_seq %in% c(rare_ASV$DNA_seq,
                         bad_length_ASV$DNA_seq,
                         n_control_ASV$DNA_seq)) %>%
  filter(sample %in% (samples_meta %>%
           filter(type == "sample") %>%
           .$sample))

samp_zero <- cdm %>%
  group_by(sample) %>%
  count(wt = count) %>%
  filter(n < 1)

asv_zero <- cdm %>%
  group_by(DNA_seq) %>%
  count(wt = count) %>%
  filter(n < 1)

cdm <- cdm %>%
  filter(!DNA_seq %in% asv_zero$DNA_seq) %>%
  filter(!sample %in% samp_zero$sample)

# CDM clustered ASV

cdm_OTU97 <- cdm %>%
  distinct(DNA_seq) %>%
  mutate(len = nchar(DNA_seq)) %>%
  group_by(len) %>%
  mutate(OTU = bioseq::seq_cluster(DNA_seq, 0.03)) %>%
  ungroup() %>%
  mutate(OTU = paste(len, OTU, sep = "_")) %>%
  dplyr::select(-len) %>%
  left_join(cdm, .)
cdm_OTU97 <- cdm_OTU97 %>%
  group_by(sample, OTU) %>%
  summarise(count = sum(count)) %>%
  left_join(group_by(cdm_OTU97, OTU) %>%
              summarise(clustered_seq = list(unique(DNA_seq)))
  )

cdm_OTU97_expanded <- cdm_OTU97 %>%
  ungroup() %>%
  distinct(OTU, .keep_all = TRUE) %>%
  unnest(clustered_seq) %>%
  left_join(taxa_id, by = c("clustered_seq" = "DNA_seq")) %>%
  select(-sample, -count, -clustered_seq)

cdm_OTU97_expanded$lowest_level <-
  cdm_OTU97_expanded[, -1] %>%
  apply(1, function(x) max(which(!is.na(x))))

cdm_OTU97_expanded <- cdm_OTU97_expanded %>%
  group_by(OTU) %>%
  filter(!duplicated(lowest_level)) %>%
  filter(lowest_level == max(lowest_level)) %>%
  mutate(lowest_level = ifelse(lowest_level > 0, lowest_level, 0))

cdm_OTU97_expanded$lowest_taxa <-
  map2_chr(1:nrow(cdm_OTU97_expanded),
           cdm_OTU97_expanded$lowest_level,
           ~ unlist(cdm_OTU97_expanded[.x, .y + 1]))


### SPECIES ####
cdm_species <- cdm %>%
  left_join(taxa_id) %>%
  filter(!is.na(species)) %>%
  group_by(species, sample) %>%
  summarise(count = sum(count)) %>%
  ungroup()

samp_zero <- cdm_species %>%
  group_by(sample) %>%
  count(wt = count) %>%
  filter(n < 1)

species_zero <- cdm_species %>%
  group_by(species) %>%
  count(wt = count) %>%
  filter(n < 1)

cdm_species <- cdm_species %>%
  filter(!species %in% species_zero$species) %>%
  filter(!sample %in% samp_zero$sample)

rar_n <- 10^4
sel_samp <- count(cdm, sample, wt = count) %>%
  filter(n > rar_n)

cdm_rar <- cdm %>%
  filter(sample %in% sel_samp$sample) %>%
  spread_cdm(sample, DNA_seq, count) %>%
  vegan::rrarefy(rar_n) %>%
  tidy_cdm(row.name = "sample", key.name = "DNA_seq", value.name = "count") %>%
  mutate(DNA_seq = as_dna(DNA_seq))


# Subsets of sites for Rhine catchment
subset_sites_rhein <-
  sites_meta %>%
  filter(main_catchment %in% c("Rhein", "Aare", "Reuss", "Limmat"))


rm(asv_zero, samp_zero, rare_ASV, bad_length_ASV, n_control_ASV, species_zero, sel_samp)
