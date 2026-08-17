library(dplyr)
library(readr)
library(glue)
library(argparser)


## Looks up one row's fine-mapping input file URIs (genotype dosages + index,
## tensorQTL permutation results, phenotype BED) from a previously-downloaded
## Terra QTL data table (see download_avtable.R) -- reading a local file here
## instead of re-querying Terra's AVTable API once per feature in a scatter.
## Column naming differs slightly between QTL types in practice (e.g. pQTL
## tables use "SubsetDosages" where eQTL/sQTL use "GenotypeDosages"), so each
## field tries a short list of known aliases.

argv <- arg_parser("Look up fine-mapping input file URIs for one feature from a downloaded QTL data table") %>%
  add_argument("--avtable-file", help = "Local TSV of the full QTL data table (see download_avtable.R)") %>%
  add_argument("--qtl-avtable", help = "Terra data table (entity type) name, e.g. eQTL_COMB_9k -- used to derive the id column name") %>%
  add_argument("--gene-name", help = "Row ID to look up (matches the <table>_id column)") %>%
  parse_args()

qtl_table <- read_tsv(argv$avtable_file, show_col_types = FALSE)

id_col <- paste0(argv$qtl_avtable, "_id")
if (!(id_col %in% colnames(qtl_table))) {
  stop(glue("Expected an id column named '{id_col}' in table '{argv$qtl_avtable}', found: {paste(colnames(qtl_table), collapse=', ')}"))
}

row <- qtl_table %>% filter(.data[[id_col]] == argv$gene_name)
if (nrow(row) != 1) {
  stop(glue("Expected exactly 1 row in '{argv$qtl_avtable}' for {id_col} == '{argv$gene_name}', found {nrow(row)}"))
}

## Not every field is present in every QTL table (e.g. some tables source
## tensorQTL permutations from a separate file rather than a table column) --
## a missing column just yields an empty string rather than aborting, since
## the caller may not even use that particular output.
pick_col <- function(candidates) {
  present <- candidates[candidates %in% colnames(row)]
  if (length(present) == 0) {
    warning(glue("None of the expected columns ({paste(candidates, collapse=', ')}) found in '{argv$qtl_avtable}' -- writing an empty value"))
    return("")
  }
  row[[present[1]]]
}

writeLines(pick_col(c("GenotypeDosages", "SubsetDosages")), "genotype_dosages.txt")
writeLines(pick_col(c("GenotypeDosagesIndex", "SubsetDosagesIndex")), "genotype_dosages_index.txt")
writeLines(pick_col(c("SubsetPermutationPvals")), "tensorqtl_permutations.txt")
writeLines(pick_col(c("SubsetBed", "SubsetPhenotypeBed")), "phenotype_bed.txt")
