library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(arrow)
library(tibble)
library(glue)
library(GenomicRanges)
library(argparser)


## Nominal cis-association test between SVA repeat-length/genotype variants
## and a molecular QTL feature set (expression, splicing, or protein), for
## individuals within `--window` bp of each other. Generic across QTL types --
## only --feature-file and --covariate-file differ between eQTL/sQTL/pQTL
## calls; --sva-file (a ProcessMEISVAData per-individual matrix) is shared.

argv <- arg_parser("Run nominal SVA repeat-length/genotype association against a QTL feature set") %>%
  add_argument("--feature-file", help = "QTL feature BED (.bed.gz) or Parquet (.parquet) file, first 4 columns chrom/start/end/feature-id") %>%
  add_argument("--sva-file", help = "SVA per-individual matrix BED (e.g. ProcessMEISVAData's deviance matrix), first 5 columns chrom/start/end/sva_id/var_id") %>%
  add_argument("--covariate-file", help = "Optional QTL covariates TSV, one row per covariate, columns are individual IDs plus an 'ID' column") %>%
  add_argument("--window", help = "cis window (bp) added on either side of each SVA locus when intersecting with features", default = 1000000) %>%
  add_argument("--prefix", help = "Output file prefix", default = "nominal_assoc") %>%
  parse_args()


## ---- Load inputs ------------------------------------------------------

if (endsWith(argv$feature_file, ".bed.gz")) {
  feature_df <- read_delim(argv$feature_file)
} else if (endsWith(argv$feature_file, ".parquet")) {
  feature_df <- read_parquet(argv$feature_file)
} else {
  stop(glue("--feature-file must end in .bed.gz or .parquet, got: {argv$feature_file}"))
}
colnames(feature_df)[1:4] <- c("chrom", "start", "end", "feature")

sva_df <- read_delim(argv$sva_file)
colnames(sva_df)[1] <- "chrom"

covariate_mat <- NULL
if (!is.na(argv$covariate_file)) {
  covariate_mat <- read_delim(argv$covariate_file) %>%
    column_to_rownames("ID") %>%
    as.matrix() %>%
    t()
}


## ---- Intersect SVA loci with features within the cis window -----------

pairfeature_intersect <- function(df1, df2, window = 1000000, id_col1 = "id1", id_col2 = "id2") {
  stopifnot(id_col1 != id_col2)
  gr1 <- GRanges(seqnames = df1$chrom,
                 ranges = IRanges(start = df1$start - window, end = df1$end + window),
                 id = df1[[id_col1]])

  gr2 <- GRanges(seqnames = df2$chrom,
                 ranges = IRanges(start = df2$start, end = df2$end),
                 id = df2[[id_col2]])

  # find overlaps
  hits <- findOverlaps(gr1, gr2)

  # extract IDs that intersect
  ids1 <- mcols(gr1)$id[queryHits(hits)]
  ids2 <- mcols(gr2)$id[subjectHits(hits)]

  df <- tibble(id1 = ids1, id2 = ids2)
  colnames(df) <- c(id_col1, id_col2)
  df
}

sva_feature_pairs <- pairfeature_intersect(
  feature_df %>% select(1:4),
  sva_df %>% select(1:4) %>% distinct(),
  window = argv$window,
  id_col1 = "feature", id_col2 = "sva_id"
)

feature_df <- feature_df %>%
  filter(feature %in% sva_feature_pairs$feature)


## ---- Build matrices for regression -------------------------------------

sva_mat <- sva_df %>%
  select(-c(1:4)) %>%
  column_to_rownames("var_id") %>%
  as.matrix()

feature_mat <- feature_df %>%
  select(-c(chrom, start, end)) %>%
  column_to_rownames("feature") %>%
  as.matrix()

intersect_indivs <- intersect(colnames(sva_mat), colnames(feature_mat))
if (!is.null(covariate_mat)) {
  intersect_indivs <- intersect(intersect_indivs, rownames(covariate_mat))
}


## ---- Run one regression per SVA x feature pair --------------------------
## `feature ~ <all SVA variant-type columns> + <all covariates>`, restricted
## to individuals with a non-missing call at every SVA variant-type column.

models <- vector("list", length = nrow(sva_feature_pairs))
for (i in seq_len(nrow(sva_feature_pairs))) {
  if (i %% 100 == 0) {
    cat(glue("\r{i}"))
  }
  if (i %% 500 == 0) {
    gc()
  }
  sva_id <- sva_feature_pairs$sva_id[i]
  feature_id <- sva_feature_pairs$feature[i]

  feature_vals <- feature_mat[feature_id, intersect_indivs]

  # Anchored so e.g. sva_id "SVA_1" doesn't also match "SVA_10_hexamer";
  # drop=FALSE so a single-variant-type match stays a matrix (not a vector
  # that t() would silently transpose into the wrong orientation).
  sva_rows <- sva_mat[grepl(paste0("^", sva_id, "_"), rownames(sva_mat)), intersect_indivs, drop = FALSE]

  data <- cbind(
    data.frame(feature = feature_vals),
    sva_rows %>% t()
  )

  test_cols <- colnames(data)[-1]

  no_na_rows <- rowSums(!is.na(data[test_cols])) == length(test_cols)
  if (sum(no_na_rows) <= 5) {
    next
  }

  if (!is.null(covariate_mat)) {
    data <- cbind(data, covariate_mat[intersect_indivs, ])
  }
  model <- lm(feature ~ ., data = data)
  models[[i]] <- coefficients(summary(model))
  rm(data)
  rm(model)
}

saveRDS(models, file = glue("{argv$prefix}_nominal_assoc.Rds"))


## ---- Extract p-values and betas for the SVA variant-type coefficients ---

all_coefs <- lapply(seq_len(nrow(sva_feature_pairs)), function(i) {
  if (i %% 100 == 0) {
    cat(glue("\r{i}"))
  }
  mat <- models[[i]]
  if (is.null(mat)) {
    return(NULL)
  }
  as.data.frame(mat) %>%
    rownames_to_column("var_id") %>%
    as_tibble() %>%
    filter(grepl("SVA", var_id)) %>%
    mutate(sva_id = sva_feature_pairs$sva_id[i], feature = sva_feature_pairs$feature[i])
}) %>%
  bind_rows() %>%
  select(sva_id, feature, var_id, beta = Estimate, p = `Pr(>|t|)`)

all_coefs %>%
  write_tsv(glue("{argv$prefix}_nominal_assoc_summary.txt.gz"))
