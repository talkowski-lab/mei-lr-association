library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(stringdist)
library(glue)
library(argparser)
library(Biostrings)
library(DECIPHER)

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
##   ins_methyl_hap_summ -- one row per (locus, haplotype), restricted to
##     loci with status "ok" and, for compound hets, to reads confidently
##     assigned to a genotyped allele. Each group's reads are aligned
##     together against the genotyped ALT allele (one multiple alignment
##     per group, so every read shares the same ALT-relative coordinate
##     system) and their per-base methylation calls are remapped from
##     insertion-sequence offsets to ALT positions, then summarized per ALT
##     position into a compact `label:pos:n_calls:n_methylated:mean:sd`
##     string, alongside mean sequence identity to ALT and the ALT
##     allele's CpG count.

argv <- arg_parser("Summarize per-read insertion methylation calls for one individual against their MEI genotypes") %>%
  add_argument("--indiv", help = "Individual ID to process (must match the `indiv` column of --genotype-table)") %>%
  add_argument("--methylation-table", help = "Per-individual insertion methylation TSV from ExtractInsertionMethylation") %>%
  add_argument("--genotype-table", help = "Combined SVA/LINE1/Alu genotype TSV (var_id, indiv, GT, hap1, hap2, ...) across all individuals; filtered internally to --indiv") %>%
  add_argument("--mei-reference-table", help = "Combined SVA/LINE1/Alu allele reference TSV (var_id, ALT_md5, ALT, length), one row per distinct insertion allele") %>%
  add_argument("--min-insertion-length", help = "Drop methylation calls for insertions shorter than this (bp)", default = 20) %>%
  add_argument("--match-bp", help = "Absolute bp difference below which a haplotype length/sequence match is accepted", default = 10) %>%
  add_argument("--match-perc", help = "Percent difference below which a haplotype length/sequence match is accepted", default = 3) %>%
  add_argument("--prob-threshold", help = "Probability >= threshold counts as methylated", default = 0.5) %>%
  add_argument("--correct-forward-strand-offset",
               help = "Apply the +1 post-hoc correction to forward-strand (`+`) methylation offsets, for methylation tables extracted before the pysam/htslib MM-tag decode fix in extract_insertion_methylation.py",
               default = TRUE) %>%
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

# Now calculate string distance and align sequences
ins_methyl_data_filtered <- ins_methyl_data %>%
  filter(var_id %in% (ins_methyl_data_summ %>% filter(status %in% c("ok", "homo_missing_one")) %>% pull(var_id))) %>%
  filter(!is.na(haplotype)) %>%
  left_join(
    ins_methyl_data_summ %>%
      mutate(gt_haplotype = case_when(
        hap1 & hap2 ~ "1,2",
        hap1 ~ "1",
        hap2 ~ "2",
        .default = NA
      )) %>%
      select(var_id, ALT_md5, gt_haplotype), relationship="many-to-many" # The many-to-many should only be from the compound hets
  ) %>%
  mutate(
    gt_haplotype = as.numeric(case_when(
      gt_haplotype == "1,2" ~ as.character(haplotype), # Arbitrary assignment if insertions are the same,
      .default = gt_haplotype
    ))
  ) %>%
  filter(
    (!(var_id %in% ins_compound_het)) |
      (assignments[var_id] == "match_hapn" & haplotype == gt_haplotype) |
      (assignments[var_id] == "antimatch_hapn" & haplotype != gt_haplotype)
  ) %>%
  distinct() %>%
  left_join(all_mei_df %>% select(-length))

# Parse "label:offset:prob,..." into a tibble of calls, applying the
# post-hoc +1 correction for forward-strand ("+") reads when the upstream
# extraction predates the extract_insertion_methylation.py pysam/htslib fix.
parse_methylation_string <- function(s, strand, correct_forward_strand_offset = argv$correct_forward_strand_offset) {
  if (is.na(s) || s == "") {
    return(tibble(offset = integer(), label = character(), prob = double()))
  }
  parts <- str_split_fixed(str_split(s, ",")[[1]], ":", 3)
  offset <- as.integer(parts[, 2])
  if (correct_forward_strand_offset && strand == "+") {
    offset <- offset + 1L
  }
  tibble(label = parts[, 1], offset = offset, prob = as.double(parts[, 3]))
}

build_group_alignment <- function(group_df) {
  seqs <- DNAStringSet(c(
    ALT = group_df$ALT[1],
    setNames(group_df$insertion_sequence, group_df$read_key)
  ))
  aln <- AlignSeqs(seqs, verbose = FALSE, processors = 1)
  strsplit(as.character(aln), "")  # named list of per-sequence aligned char vectors
}

# Per-column ungapped 0-based offset within one aligned sequence (NA at gaps).
column_to_seq_offset <- function(aligned_chars) {
  is_gap <- aligned_chars == "-"
  offset <- cumsum(!is_gap) - 1L
  offset[is_gap] <- NA_integer_
  offset
}

# offset (read coord) -> alt_position, both anchored to the same MSA columns.
offset_to_alt_map_msa <- function(aln_chars, read_key) {
  alt_off <- column_to_seq_offset(aln_chars[["ALT"]])
  read_off <- column_to_seq_offset(aln_chars[[read_key]])
  keep <- !is.na(alt_off) & !is.na(read_off)
  tibble(offset = read_off[keep], alt_position = alt_off[keep])
}

# Fraction of ALT-vs-read columns that agree, restricted to columns where
# neither is a gap (reused directly from the group MSA already computed
# in build_group_alignment -- no second alignment pass needed).
msa_identity <- function(aln_chars, read_key) {
  alt_chars <- aln_chars[["ALT"]]
  read_chars <- aln_chars[[read_key]]
  both_present <- alt_chars != "-" & read_chars != "-"
  mean(alt_chars[both_present] == read_chars[both_present])
}

summarize_group <- function(group_df) {
  group_df <- group_df %>% mutate(read_key = as.character(row_number()))
  aln_chars <- build_group_alignment(group_df)

  remapped <- group_df %>%
    mutate(
      calls = map2(insertion_methylation_string, strand, parse_methylation_string),
      align_map = map(read_key, ~ offset_to_alt_map_msa(aln_chars, .x)),
      calls = map2(calls, align_map, ~ inner_join(.x, .y, by = "offset"))
    ) %>%
    select(calls) %>%
    unnest(calls)

  pos_summary <- remapped %>%
    group_by(alt_position, label) %>%
    summarize(
      num_calls = n(),
      num_methylated = sum(prob >= argv$prob_threshold),
      methyl_fraction = num_methylated / num_calls,
      mean_prob = mean(prob),
      sd_prob = sd(prob),
      .groups = "drop"
    ) %>%
    arrange(alt_position, label)

  mean_methylation_unweighted <- with(pos_summary, mean(methyl_fraction))
  mean_methylation_weighted <- with(pos_summary, weighted.mean(methyl_fraction, num_calls))

  tibble(
    var_id = group_df$var_id[1],
    haplotype = group_df$haplotype[1],
    ALT_md5 = group_df$ALT_md5[1],
    gt_haplotype = group_df$gt_haplotype[1],
    n_reads = nrow(group_df),
    mean_alt_identity = round(mean(map_dbl(group_df$read_key, ~msa_identity(aln_chars, .x))) * 100, 2),
    num_methylation_sites = n_distinct(pos_summary$alt_position),
    mean_methylation_unweighted = round(mean_methylation_unweighted, 3),
    mean_methylation_weighted = round(mean_methylation_weighted, 3),
    num_cpg_alt = str_count(group_df$ALT[1], "CG"),
    alt_methylation_string = str_c(
      pos_summary$label, pos_summary$alt_position, pos_summary$num_calls, pos_summary$num_methylated,
      sprintf("%.3f", pos_summary$mean_prob), sprintf("%.4f", pos_summary$sd_prob), sep = ":"
    ) %>% str_c(collapse = ",")
  )
}

ins_methyl_hap_summ <- ins_methyl_data_filtered %>%
  group_split(var_id, haplotype) %>%
  map(summarize_group) %>%
  bind_rows()

ins_methyl_data_summ <- ins_methyl_data_summ %>%
  select(-hap1, -hap2) 


write_tsv(ins_methyl_data_summ, glue("{argv$prefix}.ins_methyl_data_summ.tsv"))
write_tsv(ins_methyl_hap_summ, glue("{argv$prefix}.ins_methyl_hap_summ.tsv"))
