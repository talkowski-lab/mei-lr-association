library(readr)
library(dplyr)
library(stringdist)
library(glue)
library(argparser)

## For one individual, summarizes per-read insertion methylation calls
## (from ExtractInsertionMethylation) against that individual's MEI
## genotype. Most loci have a single insertion allele, so each read's
## `haplotype` tag already tells you which allele it came from. Compound
## hets (2 distinct insertion alleles genotyped at the same locus) are
## disentangled by matching each read's observed insertion length (and,
## if that's ambiguous, its sequence) against the two alleles' expected
## length/sequence.
##
## Emits two tables:
##   ins_methyl_data_summ -- one row per locus: read/haplotype counts and
##     a genotype-consistency `status` (ok / mis_geno / nonMEI_ins / ...).
##   ins_methyl_hap_summ -- one row per (locus, haplotype) with methylation
##     summary stats, restricted to loci with status "ok" and, for
##     compound hets, to reads confidently assigned to a genotyped allele.

argv <- arg_parser("Summarize per-read insertion methylation calls for one individual against their MEI genotypes") %>%
  add_argument("--indiv", help = "Individual ID to process (must match the `indiv` column of --genotype-table)") %>%
  add_argument("--methylation-table", help = "Per-individual insertion methylation TSV from ExtractInsertionMethylation") %>%
  add_argument("--genotype-table", help = "Combined SVA/LINE1/Alu genotype TSV (var_id, indiv, GT, hap1, hap2, ...) across all individuals; filtered internally to --indiv") %>%
  add_argument("--mei-reference-table", help = "Combined SVA/LINE1/Alu allele reference TSV (var_id, ALT_md5, ALT, length), one row per distinct insertion allele") %>%
  add_argument("--min-insertion-length", help = "Drop methylation calls for insertions shorter than this (bp)", default = 20) %>%
  add_argument("--match-bp", help = "Absolute bp difference below which a haplotype length/sequence match is accepted", default = 10) %>%
  add_argument("--match-perc", help = "Percent difference below which a haplotype length/sequence match is accepted", default = 3) %>%
  add_argument("--prefix", help = "Output file prefix", default = "ins_methyl") %>%
  parse_args()

indiv <- argv$indiv
match_bp <- argv$match_bp
match_perc <- argv$match_perc

all_mei_df <- read_tsv(argv$mei_reference_table, show_col_types = FALSE) %>%
  select(var_id, ALT_md5, ALT, length) %>%
  distinct()

indiv_gt <- read_tsv(argv$genotype_table, show_col_types = FALSE) %>%
  filter(indiv == !!indiv)

ins_methyl_data <- read_tsv(argv$methylation_table, show_col_types = FALSE) %>%
  dplyr::rename(var_id = locus_name) %>%
  filter(insertion_length >= argv$min_insertion_length)

ins_methyl_data_summ <- ins_methyl_data %>%
  group_by(var_id) %>%
  summarize(
    n_NA_hap = sum(is.na(haplotype)),
    n_hap = n_distinct(haplotype, na.rm = TRUE),
    count = sum(!is.na(haplotype))
  ) %>%
  full_join(indiv_gt %>% select(var_id, ALT_md5, GT, hap1, hap2), by = "var_id") %>%
  group_by(var_id) %>%
  mutate(status = case_when(
    all(is.na(GT)) ~ "nonMEI_ins",
    all(is.na(count)) ~ "not_detected",
    all(n_NA_hap) > 0 & all(n_hap) == 0 ~ "all_na_hap",
    all(n_hap != sum(hap1 + hap2)) ~ "mis_geno",
    all(n_hap == sum(hap1 + hap2)) ~ "ok",
    .default = NA
  ))

ins_compound_het <- ins_methyl_data_summ %>%
  filter(status == "ok") %>%
  group_by(var_id) %>%
  filter(n() > 1) %>%
  pull(var_id) %>%
  unique()

# Disentangle compound hets: decide, for each of the locus's two genotyped
# alleles, which haplotype tag its reads carry.
assign_compound_het <- function(target_var_id) {
  methyl_data_subset <- ins_methyl_data %>%
    filter(var_id == !!target_var_id)

  gt_subset <- indiv_gt %>%
    filter(var_id == !!target_var_id) %>%
    left_join(all_mei_df, by = c("var_id", "ALT_md5"))

  ## First attempt by length
  hap1_l <- gt_subset %>% filter(hap1) %>% pull(length)
  hap2_l <- gt_subset %>% filter(hap2) %>% pull(length)

  length_diff_summ <- methyl_data_subset %>%
    mutate(
      length_diff1 = abs(insertion_length - hap1_l),
      length_diff2 = abs(insertion_length - hap2_l)
    ) %>%
    group_by(haplotype) %>%
    summarize(
      avg_diff_1 = mean(length_diff1),
      avg_diff_2 = mean(length_diff2)
    ) %>%
    mutate(
      avg_perc_diff_1 = avg_diff_1 / hap1_l * 100,
      avg_perc_diff_2 = avg_diff_2 / hap2_l * 100
    ) %>%
    arrange(haplotype)

  ### Scenario 1 (1=1, 2=2)
  s1 <- with(
    length_diff_summ,
    (
      ((avg_diff_1[1] < match_bp) | (avg_perc_diff_1[1] < match_perc)) &
      ((avg_diff_2[2] < match_bp) | (avg_perc_diff_2[2] < match_perc))
    )
  )
  ### Scenario 2 (1=2, 2=1)
  s2 <- with(
    length_diff_summ,
    (
      ((avg_diff_1[2] < match_bp) | (avg_perc_diff_1[2] < match_perc)) &
      ((avg_diff_2[1] < match_bp) | (avg_perc_diff_2[1] < match_perc))
    )
  )
  if (s1 & !s2) {
    return("match_hapn")
  } else if (!s1 & s2) {
    return("antimatch_hapn")
  } else if (!s1 & !s2) {
    return("no_match")
  }
  # If both plausible, disentangle with sequence diff

  hap1_seq <- gt_subset %>% filter(hap1) %>% pull(ALT)
  hap2_seq <- gt_subset %>% filter(hap2) %>% pull(ALT)

  seq_diff_summ <- methyl_data_subset %>%
    mutate(
      seq_diff1 = stringdist(insertion_sequence, hap1_seq),
      seq_diff2 = stringdist(insertion_sequence, hap2_seq)
    ) %>%
    group_by(haplotype) %>%
    summarize(
      avg_diff_1 = mean(seq_diff1),
      avg_diff_2 = mean(seq_diff2)
    ) %>%
    mutate(
      avg_perc_diff_1 = avg_diff_1 / hap1_l * 100,
      avg_perc_diff_2 = avg_diff_2 / hap2_l * 100
    ) %>%
    arrange(haplotype)

  ### Scenario 1 (1=1, 2=2)
  s1 <- with(
    seq_diff_summ,
    (
      ((avg_diff_1[1] < match_bp) | (avg_perc_diff_1[1] < match_perc)) &
      ((avg_diff_2[2] < match_bp) | (avg_perc_diff_2[2] < match_perc))
    )
  )
  ### Scenario 2 (1=2, 2=1)
  s2 <- with(
    seq_diff_summ,
    (
      ((avg_diff_1[2] < match_bp) | (avg_perc_diff_1[2] < match_perc)) &
      ((avg_diff_2[1] < match_bp) | (avg_perc_diff_2[1] < match_perc))
    )
  )
  if (s1 & !s2) {
    return("match_hapn")
  } else if (!s1 & s2) {
    return("antimatch_hapn")
  } else if (!s1 & !s2) {
    return("no_match")
  } else {
    return("ambig")
  }
}

# vapply (not sapply) so the result is always a plain character vector --
# even when ins_compound_het is empty -- and indexing it by a var_id that
# isn't a compound het correctly falls through to NA rather than erroring.
assignments <- vapply(ins_compound_het, assign_compound_het, character(1))

ins_methyl_data_summ <- ins_methyl_data_summ %>%
  mutate(status = case_when(
    assignments[var_id] == "no_match" ~ "comphet_no_match",
    assignments[var_id] == "ambig" ~ "comphet_ambig",
    .default = status
  ))

ins_methyl_hap_summ <- ins_methyl_data %>%
  filter(var_id %in% (ins_methyl_data_summ %>% filter(status == "ok") %>% pull(var_id))) %>%
  filter(!is.na(haplotype)) %>%
  group_by(var_id, haplotype) %>%
  summarize(
    n_reads = n(),
    n_reads_with_methyl = sum(insertion_num_5mC_calls > 0),
    mean_num_C = mean(insertion_num_C),
    mean_calls = mean(insertion_num_5mC_calls),
    sd_calls = sd(insertion_num_5mC_calls),
    mean_methyl_perc = mean(insertion_mean_5mC_prob, na.rm = TRUE),
    sd_methyl_perc = sd(insertion_mean_5mC_prob, na.rm = TRUE),
    .groups = "drop"
  )

ins_methyl_hap_summ <- ins_methyl_hap_summ %>%
  left_join(
    ins_methyl_data_summ %>%
      mutate(gt_haplotype = case_when(
        hap1 & hap2 ~ "1,2",
        hap1 ~ "1",
        hap2 ~ "2",
        .default = NA
      )) %>%
      select(var_id, ALT_md5, gt_haplotype),
    by = "var_id", relationship = "many-to-many" # The many-to-many should only be from the compound hets
  ) %>%
  mutate(
    gt_haplotype = as.numeric(case_when(
      gt_haplotype == "1,2" ~ as.character(haplotype), # Arbitrary assignment if insertions are the same
      .default = gt_haplotype
    ))
  ) %>%
  filter(
    (!(var_id %in% ins_compound_het)) |
      (assignments[var_id] == "match_hapn" & haplotype == gt_haplotype) |
      (assignments[var_id] == "antimatch_hapn" & haplotype != gt_haplotype)
  )

stopifnot(nrow(
  ins_methyl_hap_summ %>%
    group_by(var_id, haplotype) %>%
    filter(n() > 1)
) == 0)

ins_methyl_data_summ <- ins_methyl_data_summ %>%
  select(-hap1, -hap2) %>%
  mutate(indiv = !!indiv) %>%
  select(indiv, everything())

ins_methyl_hap_summ <- ins_methyl_hap_summ %>%
  mutate(indiv = !!indiv) %>%
  select(indiv, everything())

write_tsv(ins_methyl_data_summ, glue("{argv$prefix}.ins_methyl_data_summ.tsv"))
write_tsv(ins_methyl_hap_summ, glue("{argv$prefix}.ins_methyl_hap_summ.tsv"))
