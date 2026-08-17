library(readr)
library(dplyr)
library(stringr)
library(glue)
library(jsonlite)
library(argparser)


## Filters a nominal_association.R summary table (sva_id, feature, var_id,
## beta, p) down to nominally significant SVA x feature pairs, then joins
## each pair against a per-individual SVA genotype/repeat-length BED (e.g. a
## ProcessMEISVAData matrix) to produce one wide-format row per (variant,
## feature) combination -- every repeat-length/genotype variant type for a
## given sva_id is joined in for each feature it's significantly associated
## with. For splice (sQTL) summaries, `feature` is a leafcutter cluster ID
## (chrom:start:end:cluster:gene.version) -- it's collapsed down to the bare
## gene ID (5th colon field, version suffix dropped) so multiple clusters
## for the same gene are merged.
##
## Also emits the significant pairs grouped by feature (parallel
## feature_names.json / sva_ids.json arrays), for driving a per-feature
## scatter downstream.

argv <- arg_parser("Subset a QTL association summary to nominally significant SVA x feature pairs, and gather their genotype/repeat-length data") %>%
  add_argument("--qtl-assoc-summary", help = "nominal_association.R summary TSV (sva_id, feature, var_id, beta, p)") %>%
  add_argument("--sva-gt-bed", help = "Per-individual SVA genotype/repeat-length BED (e.g. ProcessMEISVAData matrix)") %>%
  add_argument("--signature-type", help = "QTL signature type: 'expression', 'splice', or 'protein'") %>%
  add_argument("--p-threshold", help = "Nominal p-value cutoff for significance", default = 0.01) %>%
  add_argument("--prefix", help = "Output file prefix", default = "sig_sva") %>%
  parse_args()

assoc <- read_tsv(argv$qtl_assoc_summary, show_col_types = FALSE) %>%
  filter(p < argv$p_threshold)

if (argv$signature_type == "splice") {
  assoc <- assoc %>%
    mutate(feature = str_split_i(str_split_i(feature, ":", 5), fixed("."), 1))
}

sig_pairs <- assoc %>%
  distinct(sva_id, feature)

sig_pairs %>%
  write_tsv(glue("{argv$prefix}_sig_pairs.txt"))

grouped <- sig_pairs %>%
  group_by(feature) %>%
  summarize(sva_ids = list(sva_id), .groups = "drop")

write_json(grouped$feature, "feature_names.json")
write_json(grouped$sva_ids, "sva_ids.json")

gt_bed <- read_tsv(argv$sva_gt_bed, show_col_types = FALSE)
colnames(gt_bed)[1] <- "chrom"

gt_bed <- gt_bed %>%
  mutate(variant_class = str_remove(var_id, paste0("^", sva_id, "_")))

sig_pairs %>%
  left_join(gt_bed, by = "sva_id", relationship = "many-to-many") %>%
  select(-sva_id, -chrom, -start, -end) %>%
  select(variant_id = var_id, gene_id = feature, variant_class, everything()) %>%
  write_tsv(glue("{argv$prefix}_sig_sva_gt.txt.gz"))
