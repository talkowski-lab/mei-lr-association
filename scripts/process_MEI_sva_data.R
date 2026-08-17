library(dplyr)
library(tidyr)
library(readr)
library(glue)
library(argparser)


## Combines non-reference (nrMEI, from cluster_mei_gt.R's clustered GT table)
## and reference (rMEI/ref-SVA) SVA repeat-length data into:
##   - a single tidy per-haplotype length table across both sources
##   - per-individual length + deviance-from-mean matrices (raw and
##     longest-pure-segment/LPS-adjusted), restricted to SVAs common enough
##     to be worth genotyping (AF above --common-af-threshold)

argv <- arg_parser("Combine reference and non-reference MEI SVA length data into per-individual matrices") %>%
  add_argument("--nrMEI-gt-table", help = "Clustered non-reference MEI genotype/length table (output of cluster_mei_gt.R)") %>%
  add_argument("--nrMEI-sva-length", help = "Per-nrMEI SVA repeat length table (hexamer/VNTR), keyed by fasta header") %>%
  add_argument("--nrMEI-sva-length-lps", help = "Per-nrMEI SVA longest-pure-segment (LPS) length table, keyed by fasta header") %>%
  add_argument("--rMEI-sva-length-indiv", help = "Per-individual reference MEI SVA repeat length table") %>%
  add_argument("--rMEI-sva-length-lps-indiv", help = "Per-individual reference MEI SVA LPS length table") %>%
  add_argument("--refSVA-length-table", help = "Reference SVA repeat length table (tidy, per haplotype)") %>%
  add_argument("--refSVA-bed", help = "BED of reference SVA loci (chrom, start, end, subfamily, strand, score, sva_id)") %>%
  add_argument("--n-haplotypes-nrMEI", help = "Total haplotype count (2x diploid samples) for the nrMEI cohort, used as the allele-frequency denominator", default = 12680) %>%
  add_argument("--n-haplotypes-rMEI", help = "Total haplotype count (2x diploid samples) for the rMEI cohort, used as the allele-frequency denominator", default = 12387) %>%
  add_argument("--min-hexamer-length", help = "Hexamer repeat lengths below this are treated as absent (set to 0)", default = 10) %>%
  add_argument("--min-vntr-length", help = "VNTR repeat lengths below this are treated as absent (set to 0)", default = 40) %>%
  add_argument("--common-af-threshold", help = "Minimum allele frequency for an SVA to be included in the per-individual matrices", default = 0.001) %>%
  parse_args()


## ---- Load inputs ----------------------------------------------------------

nrMEI_gt_tidy <- read_delim(argv$nrMEI_gt_table) %>%
  mutate(sva_id = paste0("SVA_MEI_cl", sprintf("%05d", cluster_id)))

nrMEI_sva_lengths <- read_delim(argv$nrMEI_sva_length) %>%
  dplyr::rename(fastaheader = ID)

nrMEI_sva_length_lps <- read_delim(argv$nrMEI_sva_length_lps) %>%
  dplyr::rename(fastaheader = ID)

rMEI_length_by_indiv <- read_delim(argv$rMEI_sva_length_indiv) %>%
  mutate(indiv = as.character(indiv))

rMEI_length_lps_by_indiv <- read_delim(argv$rMEI_sva_length_lps_indiv) %>%
  mutate(indiv = as.character(indiv))

refSVA_length_tidy <- read_delim(argv$refSVA_length_table) %>%
  mutate(indiv = as.character(indiv))

ref_sva_bed <- read_delim(argv$refSVA_bed, col_names = c("chrom", "start", "end", "subfamily", "strand", "score", "sva_id"))


## ---- nrMEI: attach repeat lengths, split by haplotype ---------------------

nrMEI_gt_tidy <- nrMEI_gt_tidy %>%
  mutate(
    fastaheader = glue("{CHROM}:{POS}-{POS+1}:{ALT_md5}")
  ) %>%
  inner_join(nrMEI_sva_lengths %>% select(fastaheader, length_hexamer, length_VNTR)) %>%
  mutate(length_hexamer = if_else(length_hexamer < argv$min_hexamer_length, 0, length_hexamer)) %>%
  mutate(length_VNTR = if_else(length_VNTR < argv$min_vntr_length, 0, length_VNTR))

nrMEI_gt_lps_tidy <- nrMEI_gt_tidy %>%
  select(-length_hexamer, -length_VNTR) %>%
  left_join(nrMEI_sva_length_lps %>% select(-starts_with("motif")))

nrMEI_gt_byhap <- bind_rows(
  nrMEI_gt_tidy %>% filter(hap1) %>% mutate(hap = "h1") %>% select(indiv, hap, CHROM, POS, sva_id, length_hexamer, length_VNTR),
  nrMEI_gt_tidy %>% filter(hap2) %>% mutate(hap = "h2") %>% select(indiv, hap, CHROM, POS, sva_id, length_hexamer, length_VNTR)
) %>%
  arrange(sva_id, indiv)

nrMEI_gt_lps_byhap <- bind_rows(
  nrMEI_gt_lps_tidy %>% filter(hap1) %>% mutate(hap = "h1") %>% select(indiv, hap, CHROM, POS, sva_id, starts_with("length")),
  nrMEI_gt_lps_tidy %>% filter(hap2) %>% mutate(hap = "h2") %>% select(indiv, hap, CHROM, POS, sva_id, starts_with("length"))
) %>%
  arrange(sva_id, indiv)


## ---- Combined ref + non-ref tidy length table ------------------------------

MEI_length_byhap <- bind_rows(
  refSVA_length_tidy %>%
    # Filter just to indivs that pass all filters
    inner_join(
      rMEI_length_by_indiv %>%
        filter(filter_GT) %>%
        select(indiv, sva_id) %>%
        distinct()
    ),
  nrMEI_gt_byhap %>% select(indiv, hap, sva_id, length_hexamer, length_VNTR)
) %>%
  arrange(sva_id)

MEI_length_byhap %>%
  write_tsv("SVA_MEI_lengths_tidy.txt.gz")


## ---- nrMEI: per-SVA summary stats + common-cluster filter -----------------

nrMEI_summary <- nrMEI_gt_byhap %>%
  group_by(sva_id) %>%
  summarize(
    CHROM = unique(CHROM),
    start = min(POS),
    end = max(POS + 1),
    mean_hexamer = mean(length_hexamer),
    mean_VNTR = mean(length_VNTR),
    AC = n(),
    AF = AC / (2 * argv$n_haplotypes_nrMEI)
  )

nrMEI_lps_summary <- nrMEI_gt_lps_byhap %>%
  group_by(sva_id) %>%
  summarize(
    mean_lps_hexamer = mean(length_lps_hexamer),
    mean_lps_VNTR_1 = mean(length_lps_VNTR_1),
    mean_lps_VNTR_2 = mean(length_lps_VNTR_2),
    mean_lps_VNTR_3 = mean(length_lps_VNTR_3)
  )

common_clusters <- nrMEI_summary %>%
  filter(AF > argv$common_af_threshold) %>%
  pull(sva_id) %>%
  unique()


## ---- nrMEI: per-haplotype deviance from the per-SVA mean -------------------

nrMEI_gt_byhap_with_dev <- nrMEI_gt_byhap %>%
  select(-CHROM, -POS) %>%
  left_join(nrMEI_summary %>% select(sva_id, starts_with("mean"))) %>%
  mutate(
    deviance_hexamer = length_hexamer - mean_hexamer,
    deviance_VNTR = length_VNTR - mean_VNTR
  ) %>%
  select(-starts_with("mean"))

nrMEI_gt_lps_byhap_with_dev <- nrMEI_gt_lps_byhap %>%
  select(-CHROM, -POS) %>%
  left_join(nrMEI_lps_summary %>% select(sva_id, starts_with("mean"))) %>%
  mutate(
    deviance_lps_hexamer = length_lps_hexamer - mean_lps_hexamer,
    deviance_lps_VNTR_1 = length_lps_VNTR_1 - mean_lps_VNTR_1,
    deviance_lps_VNTR_2 = length_lps_VNTR_2 - mean_lps_VNTR_2,
    deviance_lps_VNTR_3 = length_lps_VNTR_3 - mean_lps_VNTR_3
  ) %>%
  select(-starts_with("mean"))


## ---- nrMEI: totals per individual x SVA ------------------------------------

nrMEI_length_by_indiv <- nrMEI_gt_byhap_with_dev %>%
  group_by(indiv, sva_id) %>%
  summarize(
    n_GT_INS = n(),
    tot_length_hexamer = sum(length_hexamer), tot_length_VNTR = sum(length_VNTR),
    tot_deviance_hexamer = sum(deviance_hexamer), tot_deviance_VNTR = sum(deviance_VNTR)
  )

nrMEI_length_lps_by_indiv <- nrMEI_gt_lps_byhap_with_dev %>%
  group_by(indiv, sva_id) %>%
  summarize(
    n_GT_INS = n(),
    tot_length_lps_hexamer = sum(length_lps_hexamer), tot_deviance_lps_hexamer = sum(deviance_lps_hexamer),
    tot_length_lps_VNTR_1 = sum(length_lps_VNTR_1), tot_deviance_lps_VNTR_1 = sum(deviance_lps_VNTR_1),
    tot_length_lps_VNTR_2 = sum(length_lps_VNTR_2), tot_deviance_lps_VNTR_2 = sum(deviance_lps_VNTR_2),
    tot_length_lps_VNTR_3 = sum(length_lps_VNTR_3), tot_deviance_lps_VNTR_3 = sum(deviance_lps_VNTR_3)
  )


## ---- Generate matrices -----------------------------------------------------
## Each per-individual matrix is a BED-like table (#chrom, start, end, sva_id,
## var_id, <one column per individual>), built by pivoting a repeat-type +
## individual combo wide. pivot_wider introduces NA for individual x sva_id
## combos with no rows.
##   - nrMEI: absence means "no repeat expansion observed" for an individual
##     known to carry the insertion, so NA is filled with 0.
##   - rMEI: these SVAs are directly genotyped (incl. deletions), so a missing
##     combo means the individual was truly not genotyped -- NA is left as-is.

## -- nrMEI length + deviance --

nrMEI_length_indiv_mat <- nrMEI_length_by_indiv %>%
  filter(sva_id %in% common_clusters) %>%
  ungroup() %>%
  select(-starts_with("tot_deviance")) %>%
  pivot_longer(c(n_GT_INS, tot_length_hexamer, tot_length_VNTR), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_length_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    nrMEI_summary %>% select(CHROM, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = CHROM, start, end, sva_id, var_id, everything()) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

nrMEI_deviance_indiv_mat <- nrMEI_length_by_indiv %>%
  filter(sva_id %in% common_clusters) %>%
  ungroup() %>%
  select(-starts_with("tot_length")) %>%
  pivot_longer(c(n_GT_INS, tot_deviance_hexamer, tot_deviance_VNTR), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_deviance_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    nrMEI_summary %>% select(CHROM, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = CHROM, start, end, sva_id, var_id, everything()) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))


## -- rMEI length + deviance --

common_rMEI <- rMEI_length_by_indiv %>%
  distinct() %>%
  group_by(sva_id) %>%
  summarize(AF = sum(n_GT_INS) / (2 * argv$n_haplotypes_rMEI)) %>%
  filter(AF > argv$common_af_threshold) %>%
  pull(sva_id)

rMEI_length_indiv_mat <- rMEI_length_by_indiv %>%
  distinct() %>%
  filter(filter_GT, sva_id %in% common_rMEI) %>%
  select(indiv, sva_id, n_GT_INS, tot_length_hexamer, tot_length_VNTR) %>%
  pivot_longer(c(n_GT_INS, tot_length_hexamer, tot_length_VNTR), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_length_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    ref_sva_bed %>% select(chrom, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = chrom, start, end, sva_id, var_id, everything())

rMEI_deviance_indiv_mat <- rMEI_length_by_indiv %>%
  distinct() %>%
  filter(filter_GT, sva_id %in% common_rMEI) %>%
  select(indiv, sva_id, n_GT_INS, tot_deviance_hexamer, tot_deviance_VNTR) %>%
  pivot_longer(c(n_GT_INS, tot_deviance_hexamer, tot_deviance_VNTR), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_deviance_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    ref_sva_bed %>% select(chrom, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = chrom, start, end, sva_id, var_id, everything())

MEI_length_indiv_mat <- bind_rows(
  rMEI_length_indiv_mat,
  nrMEI_length_indiv_mat
) %>%
  arrange(`#chrom`, start, end)

MEI_deviance_indiv_mat <- bind_rows(
  rMEI_deviance_indiv_mat,
  nrMEI_deviance_indiv_mat
) %>%
  arrange(`#chrom`, start, end)

MEI_length_indiv_mat %>%
  write_tsv("mei_SVA_length_perindiv_matrix.bed.gz")

MEI_deviance_indiv_mat %>%
  write_tsv("mei_SVA_deviance_perindiv_matrix.bed.gz")


## ---- Same as above, but for LPS (longest-pure-segment) lengths ------------

## -- nrMEI LPS length + deviance --

nrMEI_length_lps_indiv_mat <- nrMEI_length_lps_by_indiv %>%
  filter(sva_id %in% common_clusters) %>%
  ungroup() %>%
  select(indiv, sva_id, n_GT_INS, starts_with("tot_length")) %>%
  pivot_longer(c(n_GT_INS, starts_with("tot_length")), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_length_lps_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    nrMEI_summary %>% select(CHROM, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = CHROM, start, end, sva_id, var_id, everything()) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

nrMEI_deviance_lps_indiv_mat <- nrMEI_length_lps_by_indiv %>%
  filter(sva_id %in% common_clusters) %>%
  ungroup() %>%
  select(indiv, sva_id, n_GT_INS, starts_with("tot_deviance")) %>%
  pivot_longer(c(n_GT_INS, starts_with("tot_deviance")), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_deviance_lps_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    nrMEI_summary %>% select(CHROM, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = CHROM, start, end, sva_id, var_id, everything()) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))


## -- rMEI LPS length + deviance --

rMEI_length_lps_indiv_mat <- rMEI_length_lps_by_indiv %>%
  distinct() %>%
  filter(filter_GT, sva_id %in% common_rMEI) %>%
  select(indiv, sva_id, n_GT_INS, starts_with("tot_length")) %>%
  pivot_longer(c(n_GT_INS, starts_with("tot_length")), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_length_lps_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    ref_sva_bed %>% select(chrom, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = chrom, start, end, sva_id, var_id, everything())

rMEI_deviance_lps_indiv_mat <- rMEI_length_lps_by_indiv %>%
  distinct() %>%
  filter(filter_GT, sva_id %in% common_rMEI) %>%
  select(indiv, sva_id, n_GT_INS, starts_with("tot_deviance")) %>%
  pivot_longer(c(n_GT_INS, starts_with("tot_deviance")), names_to = "repeat_type", values_to = "GT") %>%
  mutate(repeat_type = case_when(
    repeat_type == "n_GT_INS" ~ "nGT",
    startsWith(repeat_type, "tot") ~ sub("tot_deviance_lps_", "", repeat_type)
  )) %>%
  pivot_wider(names_from = "indiv", values_from = "GT") %>%
  left_join(
    ref_sva_bed %>% select(chrom, start, end, sva_id)
  ) %>%
  mutate(
    var_id = paste(sva_id, repeat_type, sep = "_")
  ) %>%
  select(-repeat_type) %>%
  select("#chrom" = chrom, start, end, sva_id, var_id, everything())

MEI_length_lps_indiv_mat <- bind_rows(
  rMEI_length_lps_indiv_mat,
  nrMEI_length_lps_indiv_mat
) %>%
  arrange(`#chrom`, start, end)

MEI_deviance_lps_indiv_mat <- bind_rows(
  rMEI_deviance_lps_indiv_mat,
  nrMEI_deviance_lps_indiv_mat
) %>%
  arrange(`#chrom`, start, end)

MEI_length_lps_indiv_mat %>%
  write_tsv("mei_SVA_length_lps_perindiv_matrix.bed.gz")

MEI_deviance_lps_indiv_mat %>%
  write_tsv("mei_SVA_deviance_lps_perindiv_matrix.bed.gz")
