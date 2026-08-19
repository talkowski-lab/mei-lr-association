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
## (chrom:start:end:cluster:gene.version); a bare gene ID (5th colon field,
## version suffix dropped) is derived for AVTable lookups and grouping, but
## since fine-mapping downstream expects exactly one feature id per gene,
## each gene's leafcutter feature is canonicalized to a single cluster
## (first, alphabetically, for determinism) -- a gene with multiple
## significantly-associated clusters ends up represented by just one of
## them, everywhere that gene appears.
##
## Also emits the significant pairs grouped by (bare) gene (parallel
## feature_names.json / sva_ids.json arrays), for driving a per-gene scatter
## downstream. Genes not present in the AVTable (first column) are dropped
## -- they'd otherwise fail a per-gene AVTable lookup later -- and reported
## to stderr.

argv <- arg_parser("Subset a QTL association summary to nominally significant SVA x feature pairs, and gather their genotype/repeat-length data") %>%
  add_argument("--qtl-assoc-summary", help = "nominal_association.R summary TSV (sva_id, feature, var_id, beta, p)") %>%
  add_argument("--sva-gt-bed", help = "Per-individual SVA genotype/repeat-length BED (e.g. ProcessMEISVAData matrix)") %>%
  add_argument("--avtable-file", help = "Already-downloaded QTL data table TSV (id column first) -- genes not present here are dropped") %>%
  add_argument("--signature-type", help = "QTL signature type: 'expression', 'splice', or 'protein'") %>%
  add_argument("--p-threshold", help = "Nominal p-value cutoff for significance", default = 0.01) %>%
  add_argument("--prefix", help = "Output file prefix", default = "sig_sva") %>%
  parse_args()

assoc <- read_tsv(argv$qtl_assoc_summary, show_col_types = FALSE) %>%
  filter(p < argv$p_threshold)

if (argv$signature_type == "splice") {
  assoc <- assoc %>%
    mutate(bare_gene = str_split_i(str_split_i(feature, ":", 5), fixed("."), 1))

  # Downstream fine-mapping expects exactly one leafcutter feature id per
  # gene, so canonicalize to one (first, alphabetically, for determinism)
  # cluster per gene, and use that cluster's feature id for every row of
  # that gene -- regardless of which specific cluster each row's SVA was
  # actually tested against.
  canonical_feature_by_gene <- assoc %>%
    distinct(bare_gene, feature) %>%
    arrange(bare_gene, feature) %>%
    group_by(bare_gene) %>%
    slice_head(n = 1) %>%
    ungroup()

  assoc <- assoc %>%
    select(-feature) %>%
    left_join(canonical_feature_by_gene, by = "bare_gene")
} else {
  assoc <- assoc %>%
    mutate(bare_gene = feature)
}

sig_pairs <- assoc %>%
  distinct(sva_id, feature, bare_gene)

avtable_ids <- read_tsv(argv$avtable_file, show_col_types = FALSE, col_select = 1)[[1]]

missing_genes <- setdiff(unique(sig_pairs$bare_gene), avtable_ids)
if (length(missing_genes) > 0) {
  message(glue(
    "{length(missing_genes)} significant gene(s) not present in the AVTable, dropping: ",
    "{paste(missing_genes, collapse = ', ')}"
  ))
}

sig_pairs <- sig_pairs %>%
  filter(bare_gene %in% avtable_ids)

sig_pairs %>%
  select(sva_id, feature) %>%
  write_tsv(glue("{argv$prefix}_sig_pairs.txt"))

grouped <- sig_pairs %>%
  group_by(bare_gene) %>%
  summarize(sva_ids = list(sva_id), .groups = "drop")

write_json(grouped$bare_gene, "feature_names.json")
write_json(grouped$sva_ids, "sva_ids.json")

gt_bed <- read_tsv(argv$sva_gt_bed, show_col_types = FALSE)
colnames(gt_bed)[1] <- "chrom"

gt_bed <- gt_bed %>%
  mutate(variant_class = str_remove(var_id, paste0("^", sva_id, "_")))

sig_pairs %>%
  left_join(gt_bed, by = "sva_id", relationship = "many-to-many") %>%
  select(-sva_id, -chrom, -start, -end, -bare_gene) %>%
  select(variant_id = var_id, gene_id = feature, variant_class, everything()) %>%
  write_tsv(glue("{argv$prefix}_sig_sva_gt.txt.gz"))
