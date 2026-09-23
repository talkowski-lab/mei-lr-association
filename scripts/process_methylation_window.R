library(readr)
library(dplyr)
library(glue)
library(argparser)

## Bins per-5mC-site methylation calls by signed distance from a variant and
## summarizes them per (variant, distance bin). Input is the output of
## BedtoolsIntersect run with a 5mC bedgraph as -a and a window-expanded
## variant BED as -b: each row is a 5mC site paired with a variant whose
## original (pre-window-expansion) coordinates are carried through as
## var_start/var_end, which is what distance is measured against.

argv <- arg_parser("Bin 5mC methylation calls by distance from variant and summarize per (variant, distance bin)") %>%
  add_argument("--intersect-bed", help = "BedtoolsIntersect output: 5mC bedgraph (chrom, 5mC_start, 5mC_end, 5mC_score, haplotype, coverage, est_5mC_site_count, est_unmod_site_count, discrete_5mC_score) intersected against a window-expanded variant BED (chrom, start, end, var_id, strand, var_start, var_end)") %>%
  add_argument("--indiv", help = "Individual ID to label output rows with") %>%
  add_argument("--prefix", help = "Output file prefix", default = "methylation_window") %>%
  parse_args()

dist_bins <- c(
  "ups_5-10kb", "ups_1-5kb", "ups_500-1kbp", "ups_100-500bp", "ups_0-100bp",
  "intersect",
  "downs_0-100bp", "downs_100-500bp", "downs_500-1kbp", "downs_1-5kb", "downs_5-10kb"
)

bedd <- read_tsv(
  argv$intersect_bed,
  col_names = c(
    "chrom", "5mC_start", "5mC_end", "5mC_score", "haplotype", "coverage",
    "est_5mC_site_count", "est_unmod_site_count", "discrete_5mC_score",
    "_chrom", "_start", "_end", "var_id", "strand", "var_start", "var_end"
  ),
  show_col_types = FALSE
) %>%
  select(-starts_with("_")) # duplicate chrom/start/end from the variant-BED half of the join

bedd <- bedd %>%
  mutate(
    # Signed distance from the (pre-window-expansion) variant to the 5mC site,
    # oriented by strand so downstream is positive. dist1/dist2 measure from
    # each end of the variant to each end of the site; opposite signs mean the
    # site falls inside the variant, i.e. distance 0.
    dist1 = if_else(strand == "-", var_start - `5mC_end`, `5mC_start` - var_end),
    dist2 = if_else(strand == "-", var_end - `5mC_start`, `5mC_end` - var_start),
    distance = if_else(sign(dist1) * sign(dist2) == -1, 0, sign(dist1) * pmin(abs(dist1), abs(dist2)))
  ) %>%
  mutate(dist_bin = case_when(
    distance < -5000 ~ "ups_5-10kb",
    distance < -1000 ~ "ups_1-5kb",
    distance < -500 ~ "ups_500-1kbp",
    distance < -100 ~ "ups_100-500bp",
    distance < 0 ~ "ups_0-100bp",
    distance == 0 ~ "intersect",
    distance <= 100 ~ "downs_0-100bp",
    distance <= 500 ~ "downs_100-500bp",
    distance <= 1000 ~ "downs_500-1kbp",
    distance <= 5000 ~ "downs_1-5kb",
    .default = "downs_5-10kb"
  )) %>%
  mutate(dist_bin = factor(dist_bin, levels = dist_bins, ordered = TRUE))

bedd_summ <- bedd %>%
  group_by(var_id, dist_bin) %>%
  summarize(
    n_sites = n(),
    mean_score = mean(`5mC_score`),
    sd_score = sd(`5mC_score`),
    mean_discrete_score = mean(discrete_5mC_score),
    sd_discrete_score = sd(discrete_5mC_score),
    .groups = "drop"
  ) %>%
  mutate(indiv = argv$indiv) %>%
  select(indiv, everything())

write_tsv(bedd_summ, glue("{argv$prefix}.methylation_window_summary.tsv"))
