#!/usr/bin/env Rscript
## For each SVA (ID) x repeat region (hexamer/VNTR_1/VNTR_2/VNTR_3), finds the
## longest run of a single repeated submotif -- i.e. the longest stretch of
## consecutive repeat-unit calls that all match the same sequence -- and
## reports that motif and its total length (n_repeats * motif_length).

suppressPackageStartupMessages({
  library(argparser)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(glue)
})

argv <- arg_parser("Find the longest run of a repeated submotif per SVA repeat region") %>%
  add_argument("--repeat-seq-table", "Delimited file with columns ID, region, seq (one row per repeat-unit call)") %>%
  add_argument("--prefix", "Output prefix") %>%
  parse_args()

df <- read_delim(argv$repeat_seq_table, show_col_types = FALSE) %>%
  filter(region %in% c("hex", "VNTR_1", "VNTR_2", "VNTR_3")) %>%
  mutate(match_prev = (ID == lag(ID) & region == lag(region) & seq == lag(seq))) %>%
  mutate(match_prev = replace_na(match_prev, FALSE)) %>%
  mutate(match_group = cumsum(!match_prev))

top_match <- df %>%
  group_by(ID, region, match_group) %>%
  summarize(n = n(), motif = unique(seq), .groups = "drop") %>%
  filter(!grepl("N", motif)) %>%
  mutate(motif_l = nchar(motif)) %>%
  arrange(desc(n)) %>%
  group_by(ID, region) %>%
  slice_max(order_by = tibble(n, motif_l), n = 1, with_ties = FALSE) %>%
  mutate(length_lps = n * motif_l) %>%
  select(-match_group, -n, -motif_l) %>%
  mutate(region = if_else(region == "hex", "hexamer", region))

top_match_wide <- top_match %>%
  ungroup() %>%
  pivot_wider(names_from = region, values_from = c(motif, length_lps)) %>%
  mutate(
    indiv = str_split_i(ID, "-", 1),
    hap = str_extract(ID, "asm_(h[12])", 1),
    sva_id = str_extract(ID, r"(SVA_\d{4})")
  ) %>%
  select(indiv, hap, sva_id, ends_with("hexamer"), ends_with("VNTR_1"), ends_with("VNTR_2"), ends_with("VNTR_3")) %>%
  mutate(across(starts_with("length"), ~ replace_na(.x, 0)))

write_delim(top_match_wide, glue("{argv$prefix}_lps_length.txt"))
