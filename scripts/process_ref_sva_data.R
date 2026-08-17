library(readr)
library(dplyr)
library(tidyr)
library(glue)
library(argparser)


# Parse arguments
argv <- arg_parser("Process reference SVAs to incorporate filters and export MED data") %>%
  add_argument("--ref-sva-length-table", "table with reference SVA lengths") %>%
  add_argument("--ref-sva-lps-table", "table with reference SVA longest purse segment lenghts") %>%
  add_argument("--ref-sva-flag-table", "table with SVA reference flags") %>%
  add_argument("--ref-sva-bed", "bed file with locations of reference SVAs") %>%
  add_argument("--sv-vcf", "URI of vcf with SVs (importantly deletions)") %>%
  add_argument("--l1meaid-table", "INTACT MEI output containing rMEI information") %>%
  parse_args()


# Import data
ref_sva_length_df_raw <- read_delim(argv$ref_sva_length_table) %>%
  mutate(indiv = as.character(indiv))


# Repeat-length calls below a noise floor are treated as absent (0) rather
# than as a short-but-real repeat: hexamer/VNTR segments this short are not
# reliably distinguished from assembly/alignment artifacts.
ref_sva_length_df <- ref_sva_length_df_raw %>%
  mutate(length_hexamer = if_else(length_hexamer < 10, 0, length_hexamer)) %>%
  mutate(length_VNTR = if_else(length_VNTR < 40, 0, length_VNTR))

ref_sva_length_lps_df <- read_delim(argv$ref_sva_lps_table) %>%
  select(-starts_with("motif")) %>%
  mutate(indiv = as.character(indiv)) 

# Per-haplotype QC flag -> single per-(indiv, sva_id, hap) pass/fail call.
# "Missing" (no assembly coverage, so N%/skip% couldn't even be computed) is
# folded into "Ok" here because a missing haplotype isn't itself evidence of
# a bad assembly -- it's handled downstream via the length-derived presence
# calls instead.
ref_sva_flags <- read_delim(argv$ref_sva_flag_table) %>%
  mutate(indiv = as.character(indiv)) %>%
  mutate(filter = case_when(
    is.na(N_perc_flag) & is.na(large_skip_flag) & asm_flag == "Ok" ~ "Missing",
    if_all(ends_with("flag"), ~ .x == "Ok") ~ "Ok",
    .default="Fail"
  ))

# Save filtered ref_sva_lengths
ref_sva_length_df %>% 
  left_join(ref_sva_flags) %>%
  filter(filter %in% c("Ok", "Missing")) %>%
  select(colnames(ref_sva_length_df)) %>%
  write_tsv("ref_SVA_length_filtered.txt.gz")

ref_sva_flags_h1 <- ref_sva_flags %>% filter(hap == "h1")
ref_sva_flags_h2 <- ref_sva_flags %>% filter(hap == "h2")

# h1/h2 are combined by row position (not an explicit join), so this ordering
# check has to hold or filter2 below would get silently misassigned.
stopifnot(all(ref_sva_flags_h1$indiv == ref_sva_flags_h2$indiv))
stopifnot(all(ref_sva_flags_h1$sva_id == ref_sva_flags_h2$sva_id))

# Collapse the two per-haplotype filter calls into one per-individual call.
ref_sva_flag_perindiv <- ref_sva_flags_h1 %>%
  select(indiv, sva_id, filter1=filter) %>%
  mutate(filter2=ref_sva_flags_h2$filter) %>%
  mutate(filters = case_when(
    filter1 %in% c("Ok", "Missing") & filter2 %in% c("Ok", "Missing") ~ "Ok",
    filter1 %in% c("Ok", "Missing") & filter2 == "Fail" ~ "Ok,Fail",
    filter1 == "Fail" & filter2 %in% c("Ok", "Missing") ~ "Fail,Ok",
    filter1 == "Fail" & filter2 == "Fail" ~ "Fail",
    .default=NA
  )) %>%
  select(-filter1, -filter2)

ref_sva_bed <- read_delim(argv$ref_sva_bed, col_names=c("chrom", "start", "end", "subfamily", "strand", "score", "sva_id"))

vcf_url <- argv$sv_vcf

intact_mei_med <- read_delim(argv$l1meaid_table, col_names=c("CHROM", "POS", "REF", "ALT", "ID", "ME_family", "ME_subfamily", "SVLEN"))


# Calculate deviance: for each SVA, how far a given haplotype's repeat length
# is from the mean length across "Ok"-flagged haplotypes at that same SVA
# (pooling both haplotypes into one background distribution).
ref_sva_length_summary <- ref_sva_length_df %>%
  left_join(ref_sva_flags)  %>%
  filter(filter == "Ok") %>%
  group_by(sva_id) %>%
  summarize(mean_hexamer = mean(length_hexamer), mean_VNTR = mean(length_VNTR))

ref_sva_length_df_with_dev <- ref_sva_length_df %>%
  left_join(ref_sva_length_summary) %>%
  mutate(deviance_hexamer = length_hexamer - mean_hexamer, deviance_VNTR = length_VNTR - mean_VNTR) %>%
  select(-mean_hexamer, -mean_VNTR)


ref_sva_length_lps_summary <- ref_sva_length_lps_df %>%
  left_join(ref_sva_flags)  %>%
  filter(filter == "Ok") %>%
  group_by(sva_id) %>%
  summarize(
    mean_lps_hexamer = mean(length_lps_hexamer),
    mean_lps_VNTR_1 = mean(length_lps_VNTR_1),
    mean_lps_VNTR_2 = mean(length_lps_VNTR_2),
    mean_lps_VNTR_3 = mean(length_lps_VNTR_3)
  )

ref_sva_length_lps_df_with_dev <- ref_sva_length_lps_df %>%
  left_join(ref_sva_length_lps_summary) %>%
  mutate(
    deviance_lps_hexamer = length_lps_hexamer - mean_lps_hexamer, 
    deviance_lps_VNTR_1 = length_lps_VNTR_1 - mean_lps_VNTR_1, 
    deviance_lps_VNTR_2 = length_lps_VNTR_2 - mean_lps_VNTR_2, 
    deviance_lps_VNTR_3 = length_lps_VNTR_3 - mean_lps_VNTR_3, 
  ) %>%
  select(-starts_with("mean"))

# Get MED (mobile element deletion) genotypes.
#
# The graph reference already contains the reference-SVA insertion, so an
# individual's assembly lacking that repeat shows up in the SV VCF as a
# deletion (ALT) genotype rather than as an insertion. To pull the right SV
# VCF record for each ref SVA, match L1-MEAID's SVA calls to ref_sva_bed by
# proximity (within 1kb) and pick, per sva_id, whichever candidate's SVLEN is
# closest to the ref SVA's own footprint length.
sva_med <- intact_mei_med %>%
  filter(ME_family == "SVA") %>%
  left_join(ref_sva_bed, by = c("CHROM"="chrom"), relationship="many-to-many") %>%
  filter(POS > start - 1000, POS < end + 1000) %>%
  mutate(diff_len = abs(SVLEN - (end - start))) %>%
  group_by(sva_id) %>%
  filter(diff_len == min(diff_len))

med_regions <- paste(paste(sva_med$CHROM, sva_med$POS, sep=":"), collapse=",")


token <- system(
  "gcloud auth application-default print-access-token",
  intern = TRUE
)
Sys.setenv(GCS_OAUTH_TOKEN = token)

# Requires a libcurl-enabled bcftools on PATH (see envs/Dockerfile.r_analysis)
# so it can read gs:// VCF URLs directly.
command <- glue("bcftools query -r {med_regions} -H -f '[%CHROM\t%POS\t%ID\t%SAMPLE\t%GT\n]' {vcf_url} > sva_del_gt.txt")
system(command)

med_gt_tidy <- read_delim("sva_del_gt.txt")

colnames(med_gt_tidy) <- c("CHROM", "POS", "ID", "indiv", "GT")
med_gt_tidy <- med_gt_tidy %>%
  mutate(nGT = case_when(
    GT %in% c("1|1", "1/1") ~ 2,
    GT %in% c("1|0", "0|1", "1") ~ 1,
    GT %in% c("0|0", "0/0", "0") ~ 0,
    .default = NA
  )) %>%
  left_join(sva_med %>% select(-SVLEN, -subfamily))


# Per (indiv, sva_id), sum the assembly-derived repeat lengths across
# whichever haplotype rows are present. n_GT_INS is just the haplotype-row
# count (<=2); it gets replaced below by n_pres_INS, the count of haplotypes
# where the repeat actually looks present (length > 0), once both are on hand.
med_ref_sva_length_perindiv <- ref_sva_length_df_with_dev %>%
  filter(sva_id %in% med_gt_tidy$sva_id) %>%
  group_by(indiv, sva_id) %>%
  summarize(
    n_GT_INS = n(),
    n_pres_INS = sum(length_hexamer > 0 | length_VNTR > 0),
    tot_length_hexamer = sum(length_hexamer),
    tot_length_VNTR = sum(length_VNTR),
    tot_deviance_hexamer = sum(deviance_hexamer),
    tot_deviance_VNTR = sum(deviance_VNTR),
    .groups="drop"
  )

# Cross-check the two independent haplotype-resolved data sources against
# each other: for autosomes, one haplotype carrying the ref SVA (assembly
# presence) plus the other carrying the deletion (VCF ALT genotype) should
# always sum to 2 copies total; chrX allows 1 (hemizygous) as well. filter_GT
# flags sites where that doesn't hold, i.e. the assembly and the VCF disagree.
med_gt_ref_sva_length_perindiv <- med_gt_tidy %>%
  filter(!is.na(sva_id)) %>%
  select(indiv, sva_id, CHROM, n_GT_DEL=nGT) %>%
  left_join(med_ref_sva_length_perindiv) %>%
  left_join(ref_sva_flag_perindiv) %>%
  mutate(n_GT_INS = replace_na(n_GT_INS, 0)) %>%
  mutate(n_pres_INS = replace_na(n_pres_INS, 0)) %>%
  filter(!is.na(filters)) %>%
  mutate(filter_GT = if_else(
    CHROM == "chrX", (n_pres_INS + n_GT_DEL) %in% c(1,2), n_pres_INS + n_GT_DEL == 2
  )) %>%
  filter(filters == "Ok") %>%
  select(-n_GT_INS) %>%
  dplyr::rename(n_GT_INS=n_pres_INS) %>%
  mutate(across(starts_with("tot"), ~ replace_na(.x, 0))) %>%
  arrange(indiv, sva_id)

med_gt_ref_sva_length_perindiv %>%
  distinct() %>%
  write_tsv("rMEI_SVA_totlength_perindiv_filtered.txt.gz")

med_ref_sva_length_lps_perindiv <- ref_sva_length_lps_df_with_dev %>% 
  filter(sva_id %in% med_gt_tidy$sva_id) %>%
  group_by(indiv, sva_id) %>%
  summarize(
    tot_length_lps_hexamer = sum(length_lps_hexamer),
    tot_length_lps_VNTR_1 = sum(length_lps_VNTR_1),
    tot_length_lps_VNTR_2 = sum(length_lps_VNTR_2),
    tot_length_lps_VNTR_3 = sum(length_lps_VNTR_3),
    tot_deviance_lps_hexamer = sum(deviance_lps_hexamer),
    tot_deviance_lps_VNTR_1 = sum(deviance_lps_VNTR_1),
    tot_deviance_lps_VNTR_2 = sum(deviance_lps_VNTR_2),
    tot_deviance_lps_VNTR_3 = sum(deviance_lps_VNTR_3),
    .groups="drop"
  )

med_gt_ref_sva_length_lps_perindiv <- med_gt_ref_sva_length_perindiv %>%
  select(-starts_with("tot")) %>%
  left_join(med_ref_sva_length_lps_perindiv) %>% 
  mutate(across(starts_with("tot"), ~ replace_na(.x, 0))) 

med_gt_ref_sva_length_lps_perindiv %>%
  distinct() %>%
  write_tsv("rMEI_SVA_totlength_lps_perindiv_filtered.txt.gz")




# Build a per-individual x per-SVA total-length matrix (both haplotypes
# summed), independent of the MED/VCF genotype data above -- this covers all
# "Ok"-flagged SVAs, not just ones with a matched deletion genotype.
hap1_length <- ref_sva_length_df_with_dev %>% filter(hap=="h1") %>% select(-hap)
hap2_length <- ref_sva_length_df_with_dev %>% filter(hap=="h2") %>% select(-hap)

ref_sva_perindiv <- hap1_length %>%
  full_join(hap2_length, by=c("indiv", "sva_id"), suffix=c("_h1", "_h2")) %>%
  mutate(
    tot_length_hexamer = length_hexamer_h1 + length_hexamer_h2,
    tot_length_VNTR = length_VNTR_h1 + length_VNTR_h2
  ) %>%
  select(indiv, sva_id, starts_with("tot")) %>%
  left_join(ref_sva_flag_perindiv)


ref_sva_perindiv_mat <- ref_sva_perindiv %>%
  filter(filters == "Ok") %>%
  select(-filters) %>%
  pivot_longer(c(tot_length_hexamer, tot_length_VNTR), names_to="repeat_type", values_to="tot_length", names_prefix="tot_length_") %>%
  pivot_wider(names_from=indiv, values_from="tot_length") %>%
  mutate(var_id = paste(sva_id, repeat_type, sep="_")) %>%
  select(-repeat_type) %>%
  left_join(ref_sva_bed %>% select(chrom, start, end, sva_id)) %>%
  select(`#chrom`=chrom, start, end, sva_id, var_id, everything())

ref_sva_perindiv_mat %>%
  write_tsv("ref_SVA_length_perindiv_matrix.bed.gz")


# Same matrix as above, but from the longest-pure-segment (lps) length table.
hap1_length_lps <- ref_sva_length_lps_df_with_dev %>% filter(hap=="h1") %>% select(-hap)
hap2_length_lps <- ref_sva_length_lps_df_with_dev %>% filter(hap=="h2") %>% select(-hap)

ref_sva_lps_perindiv <- hap1_length_lps %>%
  full_join(hap2_length_lps, by=c("indiv", "sva_id"), suffix=c("_h1", "_h2")) %>%
  mutate(
    tot_length_lps_hexamer = length_lps_hexamer_h1 + length_lps_hexamer_h2,
    tot_length_lps_VNTR_1 = length_lps_VNTR_1_h1 + length_lps_VNTR_1_h2,
    tot_length_lps_VNTR_2 = length_lps_VNTR_2_h1 + length_lps_VNTR_2_h2,
    tot_length_lps_VNTR_3 = length_lps_VNTR_3_h1 + length_lps_VNTR_3_h2
  ) %>%
  select(indiv, sva_id, starts_with("tot")) %>%
  left_join(ref_sva_flag_perindiv)


ref_sva_lps_perindiv_mat <- ref_sva_lps_perindiv %>%
  filter(filters == "Ok") %>%
  select(-filters) %>%
  pivot_longer(c(tot_length_lps_hexamer, tot_length_lps_VNTR_1, tot_length_lps_VNTR_2, tot_length_lps_VNTR_3), names_to="repeat_type", values_to="tot_length", names_prefix="tot_length_lps_") %>%
  pivot_wider(names_from=indiv, values_from="tot_length") %>%
  mutate(var_id = paste(sva_id, repeat_type, sep="_")) %>%
  select(-repeat_type) %>%
  left_join(ref_sva_bed %>% select(chrom, start, end, sva_id)) %>%
  select(`#chrom`=chrom, start, end, sva_id, var_id, everything())

ref_sva_lps_perindiv_mat %>%
  write_tsv("ref_SVA_length_lps_perindiv_matrix.bed.gz")

