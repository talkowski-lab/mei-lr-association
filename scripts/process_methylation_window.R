library(readr)
library(dplyr)
library(glue)
library(argparser)

## Bins per-CpG methylation calls by signed distance from a variant into
## fixed-width tiles and summarizes per (variant, haplotype, distance bin).

argv <- arg_parser("Bin per-CpG methylation calls into fixed-width tiles by distance from variant and summarize per (variant, haplotype, distance bin)") %>%
  add_argument("--intersect-bed", help = "BedtoolsIntersect output: a pb-CpG-tools pileup bed (chrom, cpg_start, cpg_end, modification_probability, haplotype, coverage, est_mod_site_count, est_unmod_site_count, discrete_modification_score) intersected against a window-expanded variant BED (chrom, start, end, var_id, variant_score, strand, var_start, var_end)") %>%
  add_argument("--bin-size", help = "Bin width in bp", default = 200) %>%
  add_argument("--flank-size", help = "Flanking size (bp) each side of the variant, as used to build --intersect-bed's variant BED -- bins only extend this far upstream/downstream", default = 5000) %>%
  add_argument("--prefix", help = "Output file prefix", default = "methylation_window") %>%
  parse_args()

bin_size <- argv$bin_size
flank_size <- argv$flank_size

bedd <- read_tsv(
  argv$intersect_bed,
  col_names = c(
    "chrom", "cpg_start", "cpg_end", "modification_probability", "haplotype",
    "coverage", "est_mod_site_count", "est_unmod_site_count",
    "discrete_modification_score",
    "_chrom", "_start", "_end", "var_id", "variant_score", "strand", "var_start", "var_end"
  ),
  show_col_types = FALSE
) %>%
  select(-starts_with("_"))

bedd <- bedd %>%
  mutate(
    dist1 = if_else(strand == "-", var_start - cpg_end, cpg_start - var_end),
    dist2 = if_else(strand == "-", var_end - cpg_start, cpg_end - var_start),
    distance = if_else(sign(dist1) * sign(dist2) == -1, 0, sign(dist1) * pmin(abs(dist1), abs(dist2)))
  ) %>%
  mutate(
    bin_edge = if_else(
      distance == 0, 0,
      sign(distance) * pmin(ceiling(abs(distance) / bin_size) * bin_size, flank_size)
    ),
    dist_bin = if_else(
      bin_edge == 0, "intersect",
      glue("{if_else(bin_edge < 0, 'ups', 'downs')}_{abs(bin_edge) - bin_size + 1}-{abs(bin_edge)}")
    )
  )

tile_label <- function(sign_prefix, edge) glue("{sign_prefix}_{edge - bin_size + 1}-{edge}")
dist_bins <- c(
  rev(tile_label("ups", seq(bin_size, flank_size, by = bin_size))),
  "intersect",
  tile_label("downs", seq(bin_size, flank_size, by = bin_size))
)
bedd <- bedd %>%
  mutate(dist_bin = factor(dist_bin, levels = dist_bins, ordered = TRUE))

bedd_summ <- bedd %>%
  group_by(var_id, haplotype, dist_bin) %>%
  summarize(
    num_cpgs = n(),
    avg_coverage = mean(coverage),
    mean_methylation_score = mean(modification_probability),
    .groups = "drop"
  )

write_tsv(bedd_summ, glue("{argv$prefix}.methylation_window_summary.tsv"))
