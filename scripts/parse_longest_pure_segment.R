#!/usr/bin/env Rscript
## Per SVA (ID) x repeat region (hexamer, VNTR_1/2/3, and the whole VNTR for
## everything but the longest pure segment):
##   - longest pure segment: the longest run of consecutive repeat-unit calls
##     with the same sequence, its motif and length (n_repeats * motif_length)
##   - purity: 1 - (summed edit cost of each unit against its sva_typer motif)
##     / (summed motif length); mismatch costs 1, an inserted or deleted base 0.5
##   - homogeneity: sum of squared unit-sequence frequencies (Simpson), 1 when
##     every unit is the same sequence; entropy: Shannon entropy of the same
##     frequencies in bits, 0 when every unit is the same sequence
##   - entropy_norm: entropy / log2(n_units), its maximum for a repeat of that
##     many units, so 0 = all units identical and 1 = every unit different;
##     NA for a single unit

suppressPackageStartupMessages({
  library(argparser)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(glue)
})

argv <- arg_parser("Longest pure segment, motif purity and homogeneity per SVA repeat region") %>%
  add_argument("--repeat-seq-table", "TSV with columns ID, region, seq (one row per repeat-unit call)") %>%
  add_argument("--prefix", "Output prefix") %>%
  parse_args()

# Must match the motifs sva_typer is run with (sva_typer src/sva.rs).
MOTIFS <- c(
  hexamer = "CCCTCT",
  VNTR_1 = "GCCTCTGCCCGGCCGCCCAGTCTGGGAAGTGAGGAGC",
  VNTR_2 = "GCCCGGCCAGCCGCCCCGTCCGGGAGGAGGTGGGGGGGTCAGCCCCC",
  VNTR_3 = "GCCGCCCCGACCGGGAAGTGAGGAGCCCCTCTGCCCG"
)
REGIONS <- c("hexamer", "VNTR_1", "VNTR_2", "VNTR_3", "VNTR")

df <- read_tsv(argv$repeat_seq_table, show_col_types = FALSE,
               col_types = cols(ID = "c", region = "c", seq = "c", .default = "?")) %>%
  filter(region %in% c("hex", "VNTR_1", "VNTR_2", "VNTR_3")) %>%
  mutate(region = if_else(region == "hex", "hexamer", region),
         seq = toupper(seq))

lps <- df %>%
  mutate(match_prev = replace_na(ID == lag(ID) & region == lag(region) & seq == lag(seq), FALSE),
         match_group = cumsum(!match_prev)) %>%
  group_by(ID, region, match_group) %>%
  summarize(n = n(), motif = first(seq), .groups = "drop") %>%
  filter(!grepl("N", motif)) %>%
  mutate(motif_l = nchar(motif)) %>%
  group_by(ID, region) %>%
  slice_max(order_by = tibble(n, motif_l), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(ID, region, motif, length_lps = n * motif_l)

# Units containing N are missing data, not impurity.
units <- df %>% filter(!grepl("N", seq))

unit_cost <- units %>%
  distinct(region, seq) %>%
  group_by(region) %>%
  mutate(ref_len = nchar(MOTIFS[[first(region)]]),
         cost = pmin(as.vector(adist(seq, MOTIFS[[first(region)]],
                                     costs = list(ins = 0.5, del = 0.5, sub = 1))),
                     ref_len)) %>%
  ungroup()

units <- left_join(units, unit_cost, by = c("region", "seq"))
units <- bind_rows(units, units %>% filter(region != "hexamer") %>% mutate(region = "VNTR"))

composition <- units %>%
  group_by(ID, region) %>%
  summarize(purity = 1 - sum(cost) / sum(ref_len),
            n_units = n(),
            .groups = "drop") %>%
  left_join(
    units %>%
      count(ID, region, seq) %>%
      group_by(ID, region) %>%
      mutate(p = n / sum(n)) %>%
      summarize(homogeneity = sum(p^2),
                entropy = -sum(p * log2(p)),
                .groups = "drop"),
    by = c("ID", "region")) %>%
  mutate(entropy_norm = if_else(n_units > 1, entropy / log2(n_units), NA_real_))

stats <- c("motif", "length_lps", "purity", "homogeneity", "entropy", "entropy_norm", "n_units")

out <- full_join(lps, composition, by = c("ID", "region")) %>%
  mutate(region = factor(region, levels = REGIONS)) %>%
  pivot_wider(names_from = region, values_from = any_of(stats), names_expand = TRUE,
              names_glue = "{.value}_{region}") %>%
  select(ID, any_of(as.vector(outer(stats, REGIONS, paste, sep = "_"))),
         -motif_VNTR, -length_lps_VNTR) %>%
  mutate(across(starts_with("length_lps_") | starts_with("n_units_"), ~ replace_na(.x, 0)))

write_tsv(out, glue("{argv$prefix}_lps_length.tsv"), na = "NA")
