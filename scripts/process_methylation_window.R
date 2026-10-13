library(readr)
library(dplyr)
library(glue)
library(argparser)

## Bins per-CpG methylation calls by signed distance from a variant into
## fixed 100bp tiles spanning the whole flanked window, and summarizes them
## per (variant, haplotype, distance bin). Input is the output of
## BedtoolsIntersect run with ONE pb-CpG-tools pileup bed (haplotype-
## specific, e.g. hap1 or hap2, or a combined/unphased file -- this script
## doesn't care which, since haplotype is read straight off the pbcpg bed's
## own embedded column rather than assumed from which file was passed in)
## as -a and a window-expanded variant BED (from CombineAndWindowVariants)
## as -b: each row is a CpG site paired with a variant whose original
## (pre-window-expansion) coordinates are carried through as
## var_start/var_end, which is what distance is measured against. Callers
## that need a per-haplotype breakdown run this once per haplotype file and
## concatenate the resulting summaries afterward, rather than concatenating
## the raw intersect beds beforehand.

argv <- arg_parser("Bin per-CpG methylation calls into 100bp tiles by distance from variant and summarize per (variant, haplotype, distance bin)") %>%
  add_argument("--intersect-bed", help = "BedtoolsIntersect output: a pb-CpG-tools pileup bed (chrom, cpg_start, cpg_end, modification_probability, haplotype, coverage, est_mod_site_count, est_unmod_site_count, discrete_modification_score) intersected against a window-expanded variant BED (chrom, start, end, var_id, strand, var_start, var_end)") %>%
  add_argument("--indiv", help = "Individual ID to label output rows with") %>%
  add_argument("--bin-size", help = "Bin width in bp", default = 100) %>%
  add_argument("--flank-size", help = "Flanking size (bp) each side of the variant, as used to build --intersect-bed's variant BED -- bins only extend this far upstream/downstream", default = 10000) %>%
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
    "_chrom", "_start", "_end", "var_id", "strand", "var_start", "var_end"
  ),
  show_col_types = FALSE
) %>%
  select(-starts_with("_")) # duplicate chrom from the variant-BED half of the join

bedd <- bedd %>%
  mutate(
    # Signed distance from the (pre-window-expansion) variant to the CpG
    # site, oriented by strand so downstream is positive. dist1/dist2
    # measure from each end of the variant to each end of the site;
    # opposite signs mean the site falls inside the variant, i.e. distance 0.
    dist1 = if_else(strand == "-", var_start - cpg_end, cpg_start - var_end),
    dist2 = if_else(strand == "-", var_end - cpg_start, cpg_end - var_start),
    distance = if_else(sign(dist1) * sign(dist2) == -1, 0, sign(dist1) * pmin(abs(dist1), abs(dist2)))
  ) %>%
  mutate(
    # Outer edge (bp) of this distance's bin_size-wide tile, capped at
    # +/-flank_size; 0 is reserved for sites inside the variant itself.
    bin_edge = if_else(
      distance == 0, 0,
      sign(distance) * pmin(ceiling(abs(distance) / bin_size) * bin_size, flank_size)
    ),
    dist_bin = if_else(
      bin_edge == 0, "intersect",
      glue("{if_else(bin_edge < 0, 'ups', 'downs')}_{abs(bin_edge) - bin_size + 1}-{abs(bin_edge)}")
    )
  )

# Order bins upstream -> intersect -> downstream, by increasing distance from
# the variant, for a readable sorted output (built programmatically since
# bin_size/flank_size determine how many tiles there are).
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
  ) %>%
  mutate(indiv = argv$indiv) %>%
  select(indiv, everything())

write_tsv(bedd_summ, glue("{argv$prefix}.methylation_window_summary.tsv"))
