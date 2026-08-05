library(argparser)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(Biostrings)
library(glue)


p <- arg_parser("Program to annotate ref SVAs with Flags for Quality") %>%
  add_argument("--sva-ref-bed", "BED file with SVA ref positions") %>%
  add_argument("--sva-ref-seqs", "FASTA file with extracted SVA ref seqs") %>%
  add_argument("--sva-ref-cov-table-hap1", "Coverage table for asm1") %>%
  add_argument("--sva-ref-cov-table-hap2", "Coverage table for asm2") %>%
  add_argument("--sva-ref-rep-coords", "Parsed repeat coordinates for SVA refs") %>%
  add_argument("--prefix", "Output Prefix") %>%
  parse_args()


ref_sva_bed_df <- read_delim(p$sva_ref_bed, col_names=c("chrom", "start", "end", "subfamily", "strand", "score", "sva_id"), delim="\t")
sva_ref_seqs <- readDNAStringSet(p$sva_ref_seqs)
hap1_cov_table <- read_delim(p$sva_ref_cov_table_hap1, delim="\t")
hap2_cov_table <- read_delim(p$sva_ref_cov_table_hap2, delim="\t")

rep_coords <- read_delim(p$sva_ref_rep_coords, delim="\t")


# Get the sequence info on the extracted ref sequences
meta_df <- tibble(
  seq = as.character(sva_ref_seqs),
  sva_seq_name = names(sva_ref_seqs)
) %>%
  mutate(
    sva_len=nchar(seq),
    n_count = str_count(seq, "N"),
    n_perc = n_count / sva_len * 100
  ) %>%
  mutate(
    indiv = str_split_i(sva_seq_name, "-", 1),
    hap = str_extract(sva_seq_name, "asm_(h[12])", 1),
    sva_id = str_extract(sva_seq_name, r"(SVA_\d{4})")
  ) %>%
  select(-sva_seq_name)

n_flag <- meta_df %>%
  group_by(sva_id, hap) %>%
  summarize(N_perc_flag=if_else(any(n_perc > 10), "High_N", "Ok"))

# Now get the status of the haplotypes
cov_table1 <- hap1_cov_table %>%
  full_join(ref_sva_bed_df %>% dplyr::rename(sva_start=start, sva_end=end), relationship="many-to-many") %>%
  filter(start > sva_start & end < sva_end)

cov_table2 <- hap2_cov_table %>%
  full_join(ref_sva_bed_df %>% dplyr::rename(sva_start=start, sva_end=end), relationship="many-to-many") %>%
  filter(start > sva_start & end < sva_end)

bad1 <- cov_table1 %>% 
  group_by(sva_id) %>%
  filter(any(footprint_depth != 1) | is.na(footprint_depth)) %>%
  mutate(mark = case_when(
    any(footprint_depth == 0) ~ "Missing",
    all(footprint_depth > 1) ~ "Remap",
    .default = "Misc"
  ))


bad2 <- cov_table2 %>% 
  group_by(sva_id) %>%
  filter(any(footprint_depth != 1) | is.na(footprint_depth)) %>%
  mutate(mark = case_when(
    any(footprint_depth == 0) ~ "Missing",
    all(footprint_depth > 1) ~ "Remap",
    .default = "Misc"
  ))

asm_bad_annot <- bad1 %>%
  select(sva_id, mark) %>%
  distinct() %>%
  mutate(hap="h1") %>%
  full_join(bad2 %>% select(sva_id, mark) %>% distinct() %>% mutate(hap="h2")) %>%
  arrange(sva_id, hap) %>%
  select(sva_id, hap, asm_flag=mark)

  # mutate(asm_flag = case_when(
  #   mark1 %in% c("Remap", "Misc") | mark2 %in% c("Remap", "Misc") ~ "Remap",
  #   mark1 == "Missing" | mark2 == "Missing" ~ "Missing",
  #   .default="Ok"
  # )) %>%

# Now look for large insertions (skip sections) within the repeats
rep_coord_summ <- rep_coords %>%
  group_by(ID, region) %>%
  summarize(total_length = sum(end-start), .groups="drop") 

rep_coord_flag_df <- rep_coord_summ %>%
  filter(region %in% c("VNTR_region", "hexamer_region")) %>%
  left_join(
    rep_coord_summ %>%
      filter(region %in% c("VNTR_region_skip", "hexamer_region_skip")) %>%
      mutate(region = gsub("_skip", "", region)) %>%
      dplyr::rename(skip_length=total_length)
  ) %>%
  mutate(skip_length = replace_na(skip_length, 0)) %>%
  mutate(skip_perc = skip_length / total_length * 100) %>%
  mutate(
    indiv = str_split_i(ID, "-", 1),
    hap = str_extract(ID, "asm_(h[12])", 1),
    sva_id = str_extract(ID, r"(SVA_\d{4})")
  ) %>%
  select(-ID) %>%
  group_by(sva_id, hap) %>%
  summarize(large_skip_flag = if_else(any(skip_perc > 20), "Large_Skip", "Ok"), .groups="drop")


flag_summ <- ref_sva_bed_df %>%
  select(chrom, start, end, sva_id) %>%
  cross_join(tibble(hap=c("h1", "h2"))) %>%
  left_join(n_flag) %>%
  left_join(asm_bad_annot) %>%
  left_join(rep_coord_flag_df) %>%
  mutate(across(ends_with("flag"), ~ replace_na(.x, "Ok")))

prefix <- p$prefix

flag_summ %>%
  group_by(sva_id) %>%
  filter(any(asm_flag %in% c("Remap", "Misc"))) %>%
  select(chrom, start, end, sva_id) %>%
  write_delim(glue("{prefix}_remap_regions.bed"), delim="\t", col_names=FALSE)

flag_summ %>%
  select(sva_id, hap, ends_with("flag")) %>%
  write_delim(glue("{prefix}_ref_sva_flags.txt"), delim="\t")

