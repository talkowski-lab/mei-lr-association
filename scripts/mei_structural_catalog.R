#!/usr/bin/env Rscript
## Catalogs the structural forms of non-reference Alu/L1/SVA insertions
## (Hancks & Kazazian 2016, Mobile DNA 7:9, Fig. 2 panels b-i) from the
## RepeatMasker hits L1ME-AID stores per sequence, and optionally traces the
## genomic sources of transduced sequence. Writes tables only; figures and
## tests live in notebooks/mei_structural_catalog_analysis.Rmd.
##
## Outputs (<outdir>/<prefix>.structure_*):
##   per_sequence.tsv.gz       one row per input sequence, every call
##   per_variant.tsv           rollup per variant: modal structure, majority flags, medians, concordance
##   transductions.fa          one transduced segment per variant and kind (variant|kind|family|chrom|pos)
##   transduction_sources.tsv  source tracing per segment     (only with --reference / --transduction-paf)
##   source_summary.tsv        offspring per source locus     (same)
##   transductions.tsv         one row per transduced segment: sequence, its sub-segments and their sources, flag (same)
##   transduction_segments.tsv one row per sub-segment, innermost (next to the element) first, with its alignment (same)
##
## Inputs. The catalog starts from one row per insertion sequence plus the
## RepeatMasker hits that annotate it, read from L1ME-AID output TSVs
## (github.com/Markloftus/L1ME-AID) or the INTACT_MEI `*.limeaid.v134.intact.tsv`
## (same table, FLAGS column overwritten with INTACT / INTACT_3end). Per-hit
## query and consensus coordinates exist only in `Element_Hits` (v1.4) /
## `TE_Hits` (older), stored as a Python-list repr of raw RepeatMasker lines --
## everything downstream is built on those. Columns are looked up by header
## name because their names/order drift between L1ME-AID versions. Each
## sequence is joined to the clustered allele table on chrom + POS + md5 of the
## sequence, which gives its variant and family (MEI_type), and then to the
## genotype table on that allele's CHROM + POS + ALT_md5.
##
## Everything is dplyr: the largest tables (RepeatMasker hits, genotypes, the
## reference rmsk) are a few million rows at most, which dplyr handles in
## seconds.

suppressPackageStartupMessages({
  library(argparser)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringi)
})

## ---- Arguments -----------------------------------------------------------------

ARGS <- list(
  limeaid_table = c("", "Required. L1ME-AID/INTACT_MEI TSV covering every family"),
  mei_table = c("", "Required. Clustered MEI allele table (chrom, start, ..., ALT, MEI_type, ..., cluster_id, var_id, ALT_md5); groups sequences into variants and sets their family"),
  mei_variant_col = c("var_id", "mei-table column that defines a variant"),
  genotype_table = c("", "Required. Per-individual genotypes (CHROM, POS, ALT_md5, indiv, GT, hap1, hap2, ...), keyed like --mei-table"),
  n_individuals = c("0", "Cohort size for AF (AN = 2 x this); 0 = distinct indiv in --genotype-table"),
  max_5p_missing = c("Alu=20,L1=50,SVA=50", "Consensus bp missing at the 5' end before calling 5' truncation, per family"),
  max_3p_missing = c("Alu=20,L1=50,SVA=50", "Consensus bp missing at the 3' end before calling 3' truncation, per family"),
  min_transduction = c("30", "Minimum unique bp for a transduction"),
  f1_min_mast2_bp = c("60", "MAST2 exon 1 bp in front of an SVA that make it SVA_F1 (MAST2 then counts as element, not 5' transduction)"),
  sva_hexamer_cons_bp = c("80", "SVA consensus bp covered by its leading hexamer: an SVA hit starting within it, with hexamer directly in front, is not 5' truncated"),
  max_transduction_masked = c("0.5", "Maximum repeat-masked fraction of a transduced segment"),
  min_orphan = c("100", "Minimum unique bp for an orphan-transduction candidate"),
  min_polya = c("8", "Minimum poly(A) length"),
  polya_min_frac = c("0.8", "Minimum A fraction of a poly(A) segment"),
  polya_gap = c("20", "Maximum gap between element 3' end and its first poly(A)"),
  terminal_slack = c("30", "End bases never treated as transduced (TSD / anchor base)"),
  min_inversion = c("20", "Minimum opposite-strand bp for a 5' inversion"),
  inv_slack = c("30", "Allowed overlap between inverted and sense segments"),
  min_internal_gap = c("50", "Minimum consensus gap between sense hits for an internal deletion"),
  reference = c("", "GRCh38 FASTA or minimap2 .mmi; enables transduction source tracing"),
  reference_fasta = c("", "GRCh38 FASTA for the source-flank search (default: --reference unless it is a .mmi)"),
  reference_fasta_index = c("", "Its .fai, when not beside it (built on the fly otherwise)"),
  transduction_paf = c("", "Precomputed PAF of transductions.fa vs the reference (instead of --reference)"),
  ref_rmsk = c("", "Reference repeats: UCSC rmsk.txt[.gz], RepeatMasker .out, or BED"),
  known_sources = c("", "Comma-separated BEDs of extra (e.g. non-reference) source elements"),
  genes = c("", "Gene/exon BED to annotate source loci"),
  minimap2_preset = c("sr", "minimap2 -x preset"),
  threads = c("4", "minimap2 threads"),
  secondary_ratio = c("0.05", "minimap2 -p: keep secondary hits scoring this fraction of the primary (finds short inherited pieces of concatenated transductions)"),
  min_mapq = c("20", "MAPQ for a unique source mapping"),
  min_identity = c("0.9", "Minimum identity of a source hit"),
  min_qcov = c("0.7", "Minimum query coverage of a source hit"),
  self_window = c("1000", "Hits this close to the insertion site are self hits"),
  source_window = c("2000", "Maximum distance from a source element"),
  overlap_tolerance = c("50", "Allowed overlap with a source element"),
  cluster_window = c("5000", "Window for grouping element-less source loci"),
  full_length_l1_bp = c("5900", "L1 source length counted as full-length"),
  run_ava = c("true", "All-vs-all cluster unmapped segments into shared non-reference sources"),
  min_segment = c("20", "Shortest sub-segment of a transduction kept, and shortest untraced stretch searched against source flanks"),
  flank_search = c("true", "BLAST untraced stretches against the flanks (source_window bp) of candidate source elements"),
  flank_families = c("SVA,L1", "Families whose elements' flanks are searched (reference L1s only if full-length); Alu flanks are most of the genome"),
  flank_max_evalue = c("1e-5", "Maximum blastn-short e-value of a source-flank hit"),
  repeat_max_div = c("0.15", "A traced genome hit needs min_segment bp outside RefRmsk repeats below this divergence, Alus and simple/low-complexity sequence (the flank library is masked for every repeat)"),
  outdir = c(".", "Output directory"),
  prefix = c("mei", "Output file prefix")
)
parser <- arg_parser("Catalog MEI structural forms and trace transduction sources")
for (name in names(ARGS)) {
  parser <- add_argument(parser, paste0("--", gsub("_", "-", name)), help = ARGS[[name]][2], default = ARGS[[name]][1])
}
argv <- parse_args(parser)
## argparser hands every option over as a string; "" means unset
params <- setNames(lapply(names(ARGS), function(name) {
  value <- argv[[name]]
  if (is.null(value) || (length(value) == 1 && is.na(value))) "" else as.character(value)
}), names(ARGS))
for (required in c("limeaid_table", "mei_table", "genotype_table")) {
  if (!nzchar(params[[required]])) stop("--", gsub("_", "-", required), " is required")
}

## "Alu=20,L1=50" -> c(Alu = 20L, L1 = 50L)
per_family <- function(x) {
  pairs <- stri_split_fixed(stri_split_fixed(x, ",")[[1]], "=", n = 2)
  setNames(as.integer(stri_trim_both(map_chr(pairs, 2))), stri_trim_both(map_chr(pairs, 1)))
}
opt <- list(
  max_5p_missing = per_family(params$max_5p_missing),
  max_3p_missing = per_family(params$max_3p_missing),
  min_transduction = as.integer(params$min_transduction),
  max_transduction_masked = as.numeric(params$max_transduction_masked),
  min_orphan = as.integer(params$min_orphan),
  min_polya = as.integer(params$min_polya),
  polya_min_frac = as.numeric(params$polya_min_frac),
  polya_gap = as.integer(params$polya_gap),
  terminal_slack = as.integer(params$terminal_slack),
  min_inversion = as.integer(params$min_inversion),
  inv_slack = as.integer(params$inv_slack),
  min_internal_gap = as.integer(params$min_internal_gap),
  min_segment = as.integer(params$min_segment),
  sva_hexamer_cons_bp = as.integer(params$sva_hexamer_cons_bp),
  f1_min_mast2_bp = as.integer(params$f1_min_mast2_bp),
  min_mapq = as.integer(params$min_mapq),
  min_identity = as.numeric(params$min_identity),
  min_qcov = as.numeric(params$min_qcov),
  self_window = as.integer(params$self_window),
  source_window = as.integer(params$source_window),
  overlap_tolerance = as.integer(params$overlap_tolerance),
  cluster_window = as.integer(params$cluster_window),
  full_length_l1_bp = as.integer(params$full_length_l1_bp),
  repeat_max_div = as.numeric(params$repeat_max_div),
  run_ava = as.logical(params$run_ava) %in% TRUE,
  flank_search = as.logical(params$flank_search) %in% TRUE
)

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
out_file <- function(name) file.path(params$outdir, paste0(params$prefix, ".structure_", name))

## ---- Startup checks -------------------------------------------------------------
## The genome pass leaves short and low-MAPQ pieces to the source-flank search,
## so a run without it loses many correct sources: fail before the slow part.

run_sources <- nzchar(params$reference) || nzchar(params$transduction_paf)
reference_fasta <- if (nzchar(params$reference_fasta)) params$reference_fasta else
  if (nzchar(params$reference) && !grepl("\\.mmi$", params$reference)) params$reference else ""
if (run_sources && opt$flank_search) {
  problems <- c(
    if (!nzchar(reference_fasta)) "no reference FASTA: --reference is a .mmi or absent (only --transduction-paf); pass --reference-fasta",
    if (nzchar(reference_fasta) && !file.exists(reference_fasta)) sprintf("reference FASTA %s does not exist", reference_fasta),
    if (!nzchar(params$ref_rmsk)) "no --ref-rmsk: without it the flank library and hits cannot be repeat-masked",
    if (!nzchar(Sys.which("samtools"))) "samtools not on PATH",
    if (!nzchar(Sys.which("makeblastdb"))) "makeblastdb not on PATH",
    if (!nzchar(Sys.which("blastn"))) "blastn not on PATH")
  if (length(problems)) stop("source-flank search cannot run:\n  - ", paste(problems, collapse = "\n  - "),
                             "\nFix these, or pass --flank-search false to accept that short and low-MAPQ transduction pieces stay untraced.")
} else if (run_sources) {
  message(strrep("!", 78), "\nWARNING: --flank-search false. Short (< ~50 bp) and low-MAPQ transduction pieces will stay\n",
          "untraced (ambiguous_source / unaligned); in simulation 50 bp recall fell from ~92% to ~61%.\n", strrep("!", 78))
}

## ---- L1ME-AID table: insertions and their RepeatMasker hits ----------------------
## L1ME-AID's ID is `chrom:POS;REF`, which two different ALT alleles at one site
## share, so every row gets its own key (`id`) and hits are keyed on it; the
## original ID is kept as `limeaid_id` for the allele-table join. A headerless
## table (shards concatenated without one, as AnnotateL1MEAID's
## ConcatIntactMEIShards does) is read in the L1ME-AID v1.4 column order. Older
## versions name two of these TE_Hits / FILTER_RESULTS instead.

LIMEAID_COLS <- c("ID", "Sequence", "Element_Hits", "Sequence_Length", "Element_Annotation",
                  "Element_Divergence", "Element_Proportion", "Element_Percentage", "Element_Designation",
                  "Orientation", "FLAGS", "Tail_Begins", "Tail_Type", "Tail_Length", "Tail_Seed_Hits",
                  "Unique_Element_Count", "Twin_Priming_Flag")
limeaid <- read_tsv(params$limeaid_table, col_names = FALSE, col_types = cols(.default = col_character()),
                    quote = "", na = "NA", progress = FALSE)
if (nrow(limeaid) > 0 && any(tolower(unlist(limeaid[1, ])) == "sequence")) {
  names(limeaid) <- make.unique(as.character(unlist(limeaid[1, ])))
  limeaid <- limeaid[-1, ]
} else {
  n_named <- min(ncol(limeaid), length(LIMEAID_COLS))
  names(limeaid)[seq_len(n_named)] <- LIMEAID_COLS[seq_len(n_named)]
  message(params$limeaid_table, ": no header; reading it as the L1ME-AID v1.4 column order (", ncol(limeaid), " columns)")
}
layout <- paste0("expected L1ME-AID output, either with a header or headerless in the v1.4 order: ",
                 paste(LIMEAID_COLS, collapse = ", "))
column_named <- function(...) {
  for (candidate in c(...)) {
    found <- names(limeaid)[tolower(names(limeaid)) == tolower(candidate)]
    if (length(found)) return(found[1])
  }
  NULL
}
id_col <- column_named("ID") %||% names(limeaid)[1]
seq_col <- column_named("Sequence")
hits_col <- column_named("Element_Hits", "TE_Hits")
if (is.null(seq_col) || is.null(hits_col)) stop(params$limeaid_table, ": no Sequence / Element_Hits columns; ", layout)
## Positions are only an assumption without a header, so check the content.
seq_ok <- mean(stri_detect_regex(head(limeaid[[seq_col]], 2000), "^[ACGTNacgtn]+$"))
hits_ok <- mean(stri_detect_regex(head(limeaid[[hits_col]], 2000), "^\\s*\\[") |
                  head(limeaid[[hits_col]], 2000) %in% c("", "NONE", "nan", "NA", "[]"))
if (seq_ok < 0.95 || hits_ok < 0.95) {
  stop(sprintf("%s: column '%s' is %.0f%% DNA and column '%s' is %.0f%% RepeatMasker hit lists (need >= 95%% each); %s",
               params$limeaid_table, seq_col, 100 * seq_ok, hits_col, 100 * hits_ok, layout))
}
annotation_col <- column_named("Element_Annotation")
flag_col <- column_named("FLAGS", "FILTER_RESULTS")
flags <- if (!is.null(flag_col)) limeaid[[flag_col]] else NA_character_
is_intact_flag <- stri_detect_regex(flags, "^(INTACT|soloLTR|ERV)") %in% TRUE
insertions <- tibble(
  id = paste0("row", seq_len(nrow(limeaid))),
  limeaid_id = limeaid[[id_col]],
  seq = toupper(limeaid[[seq_col]]),
  subfamily = if (!is.null(annotation_col)) limeaid[[annotation_col]] else NA_character_,
  limeaid_flags = if_else(is_intact_flag, NA_character_, flags),
  intact_flag = if_else(is_intact_flag, flags, NA_character_)
) %>%
  mutate(L = nchar(seq))

## Each Element_Hits cell is a Python-list repr of RepeatMasker .out lines.
## RepeatMasker coordinates are 1-based inclusive; on 'C'-strand hits the
## consensus columns are "(left) end begin", so they are re-ordered and
## cons_begin <= cons_end always holds.
hit_cells <- limeaid[[hits_col]]
hit_cells[is.na(hit_cells) | hit_cells %in% c("", "NONE", "nan", "NA")] <- ""
quoted <- stri_extract_all_regex(hit_cells, "'[^']*'|\"[^\"]*\"")
n_lines <- lengths(quoted)
n_lines[map_lgl(quoted, ~ length(.x) == 1 && is.na(.x[1]))] <- 0L
rm_lines <- stri_trim_both(stri_replace_all_fixed(stri_sub(unlist(quoted[n_lines > 0]), 2, -2), "*", " "))
tokens <- stri_split_regex(rm_lines, "\\s+")
has_fields <- lengths(tokens) >= 14
fields <- do.call(rbind, lapply(tokens[has_fields], `[`, 1:14))
if (is.null(fields)) fields <- matrix(character(), 0, 14)
unparen <- function(x) suppressWarnings(as.integer(stri_replace_all_regex(x, "[()]", "")))
on_plus <- fields[, 9] == "+"
rm_hits <- tibble(
  id = rep(insertions$id, n_lines)[has_fields],
  divergence = suppressWarnings(as.numeric(fields[, 2])),
  query_start = as.integer(fields[, 6]), query_end = as.integer(fields[, 7]),
  strand = if_else(on_plus, "+", "C"), repeat_name = fields[, 10], repeat_class = fields[, 11],
  col12 = unparen(fields[, 12]), col13 = unparen(fields[, 13]), col14 = unparen(fields[, 14])
) %>%
  mutate(cons_begin = if_else(on_plus, col12, col14),
         cons_end = col13,
         cons_left = if_else(on_plus, col14, col12),
         lo = pmin(cons_begin, cons_end), hi = pmax(cons_begin, cons_end),
         cons_begin = lo, cons_end = hi) %>%
  select(-col12, -col13, -col14, -lo, -hi) %>%
  filter(!is.na(query_start), !is.na(query_end), query_end >= query_start, !is.na(cons_begin)) %>%
  mutate(rm_family = case_when(stri_detect_regex(repeat_class, "^SINE/Alu", case_insensitive = TRUE) ~ "Alu",
                               stri_detect_regex(repeat_class, "^LINE/L1", case_insensitive = TRUE) ~ "L1",
                               stri_detect_regex(repeat_class, "(^|/)SVA", case_insensitive = TRUE) ~ "SVA",
                               TRUE ~ NA_character_))

## ---- MEI allele table: variant and family of each sequence -------------------------
## The `*.clustered.bed` lineage: one row per distinct ALT allele at a site,
## `cluster_id`/`var_id` grouping alleles into one insertion. Joined to
## L1ME-AID rows on chrom + POS + md5(ALT), because L1ME-AID IDs are
## `chrom:POS;REF` and its Sequence is the VCF ALT; md5 alone collides across
## sites (identical full-length alleles recur genome-wide). The table's own
## ALT_md5 (md5 of the ALT as written) is what genotype tables are keyed on.

MEI_TABLE_COLS <- c("chrom", "start", "end", "REF", "ALT", "MEI_type", "MEI_subfamily", "length",
                    "orientation", "intact_status", "cluster_id", "var_id", "ALT_md5")
mei_table <- read_tsv(params$mei_table, col_types = cols(.default = col_character()), na = "NA", progress = FALSE)
if (!all(c("chrom", "start") %in% names(mei_table))) {
  if (ncol(mei_table) != length(MEI_TABLE_COLS)) {
    stop(params$mei_table, ": no header and ", ncol(mei_table), " columns; expected a header or the ",
         length(MEI_TABLE_COLS), " columns ", paste(MEI_TABLE_COLS, collapse = ", "))
  }
  mei_table <- read_tsv(params$mei_table, col_names = MEI_TABLE_COLS, col_types = cols(.default = col_character()),
                        na = "NA", progress = FALSE)
}
variant_col <- params$mei_variant_col
if (!variant_col %in% names(mei_table)) stop("mei_table has no '", variant_col, "' column (mei_variant_col)")
if (!"MEI_type" %in% names(mei_table)) mei_table$MEI_type <- NA_character_
mei_table <- mei_table %>%
  mutate(md5 = if ("ALT" %in% names(.)) as.character(openssl::md5(toupper(ALT))) else tolower(ALT_md5),
         ALT_md5 = if ("ALT_md5" %in% names(.)) tolower(ALT_md5) else md5,
         pos = as.integer(start),
         variant_id = if (variant_col == "cluster_id") paste(MEI_type, cluster_id, sep = "_cl") else .data[[variant_col]])
allele_keys <- mei_table %>% distinct(chrom, pos, md5, variant_id, ALT_md5, MEI_type)
n_conflict <- allele_keys %>% group_by(chrom, pos, md5) %>% filter(n_distinct(variant_id) > 1) %>% n_groups()
allele_keys <- allele_keys %>% distinct(chrom, pos, md5, .keep_all = TRUE)

id_location <- stri_match_first_regex(insertions$limeaid_id, "^(.+?):(\\d+)(?:;|$)")
insertions <- insertions %>%
  mutate(chrom = id_location[, 2], id_pos = as.integer(id_location[, 3]), md5 = as.character(openssl::md5(seq)))
if (all(is.na(insertions$id_pos))) {
  stop("L1ME-AID IDs are not chrom:POS[;REF] (e.g. ", insertions$limeaid_id[1], "); cannot join mei_table")
}
## Table start is VCF POS in cluster_mei_gt.R; also try 0-based BED starts.
matches_by_offset <- map_int(c(0L, -1L), function(offset)
  nrow(inner_join(insertions %>% transmute(chrom, pos = id_pos + offset, md5), allele_keys, by = c("chrom", "pos", "md5"))))
pos_offset <- c(0L, -1L)[which.max(matches_by_offset)]
insertions <- insertions %>%
  mutate(pos = id_pos + pos_offset) %>%
  left_join(allele_keys, by = c("chrom", "pos", "md5"))
cat(sprintf("mei_table: %d allele rows, %d variants; joined with start = POS%s.\n",
            nrow(mei_table), n_distinct(mei_table$variant_id), if (pos_offset) " - 1" else ""))
if (n_conflict) cat(sprintf("mei_table: %d alleles map to more than one variant; kept the first.\n", n_conflict))
print(knitr::kable(insertions %>% summarise(l1meaid_rows = n(), matched = sum(!is.na(variant_id)),
                                            variants = n_distinct(variant_id[!is.na(variant_id)])),
                   caption = "L1ME-AID rows matched to mei_table (unmatched rows are dropped)"))
if (!any(!is.na(insertions$variant_id))) stop("No L1ME-AID row matched mei_table on chrom + POS + md5(ALT)")
insertions <- insertions %>% filter(!is.na(variant_id)) %>% select(-id_pos, -md5)

## Family from the allele table's MEI_type, else the RepeatMasker family with
## the most bp in the sequence; sequences with neither are dropped. A blank
## L1ME-AID subfamily is filled from the longest same-family hit.
dominant_family <- rm_hits %>%
  filter(!is.na(rm_family)) %>%
  mutate(seen = row_number()) %>%
  group_by(id, rm_family) %>%
  summarise(bp = sum(query_end - query_start + 1L), seen = min(seen), .groups = "drop") %>%
  arrange(id, desc(bp), seen) %>%
  distinct(id, .keep_all = TRUE)
insertions <- insertions %>%
  mutate(family = case_when(stri_detect_regex(MEI_type, "alu", case_insensitive = TRUE) ~ "Alu",
                            stri_detect_regex(MEI_type, "line|l1", case_insensitive = TRUE) ~ "L1",
                            stri_detect_regex(MEI_type, "sva", case_insensitive = TRUE) ~ "SVA",
                            TRUE ~ NA_character_))
n_from_type <- sum(!is.na(insertions$family))
insertions <- insertions %>%
  left_join(dominant_family %>% select(id, rm_family), by = "id") %>%
  mutate(family = coalesce(family, rm_family)) %>%
  select(-rm_family)
n_without <- sum(is.na(insertions$family))
message(sprintf("family: %d from MEI_type, %d from RepeatMasker hits, %d with neither (dropped).",
                n_from_type, nrow(insertions) - n_from_type - n_without, n_without))
insertions <- insertions %>% filter(!is.na(family))
longest_family_hit <- rm_hits %>%
  inner_join(insertions %>% select(id, family), by = "id") %>%
  filter(rm_family == family) %>%
  arrange(desc(query_end - query_start)) %>%
  distinct(id, .keep_all = TRUE) %>%
  select(id, longest_hit = repeat_name)
insertions <- insertions %>%
  left_join(longest_family_hit, by = "id") %>%
  mutate(subfamily = if_else(is.na(subfamily) | subfamily %in% c("", "NONE", "nan"), longest_hit, subfamily)) %>%
  select(-longest_hit)
max_5p_missing <- opt$max_5p_missing
max_3p_missing <- opt$max_3p_missing
for (fam in setdiff(unique(insertions$family), names(max_5p_missing))) max_5p_missing[fam] <- 50L
for (fam in setdiff(unique(insertions$family), names(max_3p_missing))) max_3p_missing[fam] <- 50L

## Identical (family, sequence, hits) inputs give identical calls: classify
## each distinct one once. The hit signature leaves out the query name.
hit_signature <- rm_hits %>%
  group_by(id) %>%
  summarise(sig = paste(query_start, query_end, strand, repeat_name, cons_begin, cons_end, collapse = ";"))
insertions <- insertions %>%
  left_join(hit_signature, by = "id") %>%
  mutate(sig = coalesce(sig, ""),
         key = paste(family, seq, sig, sep = "\r"),
         uid = match(key, unique(key))) %>%
  select(-sig, -key)
representatives <- insertions %>% distinct(uid, .keep_all = TRUE) %>% select(uid, id)
message(sprintf("Loaded %d sequences (%d distinct) across %s.", nrow(insertions), nrow(representatives),
                paste(sprintf("%s: %d", names(table(insertions$family)), table(insertions$family)), collapse = ", ")))

## ---- Element-sense orientation ------------------------------------------------------
## The strand carrying most family bp becomes '+': the sequence is reverse
## complemented and every hit re-expressed in that orientation. Family hits
## are the family's own (core) hits plus, for SVA, (CCCTCT)n-style hexamer hits.

SVA_HEXAMER_RE <- "^\\((CCCTCT|CCTCTC|CTCTCC|TCTCCC|CTCCCT|TCCCTC|AGAGGG|GAGGGA|AGGGAG|GGGAGA|GGAGAG|GAGAGG)\\)n$"
TAIL_REPEAT_RE <- "^\\((A|AA[A-Z]?|A{1,3}[CGT]?|T|TT[A-Z]?|T{1,3}[ACG]?)\\)n$|^[AT]-rich$"
to_classify <- insertions %>% semi_join(representatives, by = "id")
sense_hits <- rm_hits %>%
  inner_join(to_classify %>% select(id, family, L), by = "id") %>%
  mutate(is_core = !is.na(rm_family) & rm_family == family,
         is_family = is_core | (family == "SVA" & stri_detect_regex(repeat_name, SVA_HEXAMER_RE, case_insensitive = TRUE)))
majority_strand <- sense_hits %>%
  filter(is_core) %>%
  mutate(seen = row_number()) %>%
  group_by(id, strand) %>%
  summarise(bp = sum(query_end - query_start + 1), seen = min(seen), .groups = "drop") %>%
  arrange(id, desc(bp), seen) %>%
  distinct(id, .keep_all = TRUE) %>%
  transmute(id, orientation = if_else(strand == "+", "+", "-"))
to_classify <- to_classify %>%
  left_join(majority_strand, by = "id") %>%
  mutate(strand = coalesce(orientation, "+"),
         sense = if_else(strand == "-", stri_reverse(chartr("ACGTNacgtn", "TGCANtgcan", seq)), seq)) %>%
  select(-orientation)
minus_ids <- to_classify$id[to_classify$strand == "-"]
sense_hits <- sense_hits %>%
  mutate(flip = id %in% minus_ids,
         new_start = if_else(flip, L - query_end + 1L, query_start),
         query_end = if_else(flip, L - query_start + 1L, query_end),
         query_start = new_start,
         strand = if_else(flip, if_else(strand == "+", "C", "+"), strand),
         query_len = query_end - query_start + 1L) %>%
  select(-flip, -new_start)
core_hits <- sense_hits %>% filter(is_family, is_core)
sense_core <- core_hits %>% filter(strand == "+") %>% arrange(id, query_start)
inverted_core <- core_hits %>% filter(strand == "C")
hexamer_hits <- sense_hits %>% filter(is_family, !is_core)
other_hits <- sense_hits %>% filter(!is_family)
is_tail_hit <- stri_detect_regex(other_hits$repeat_name, TAIL_REPEAT_RE, case_insensitive = TRUE)
tail_hits <- other_hits[is_tail_hit & other_hits$strand == "+", ]
nontail_hits <- other_hits[!is_tail_hit, ]

## ---- Classification: element span -------------------------------------------------

calls <- to_classify %>%
  select(id, family, subfamily, L, strand, limeaid_flags, intact_flag, sense) %>%
  left_join(core_hits %>% group_by(id) %>% summarise(family_bp = sum(query_len)), by = "id") %>%
  mutate(family_bp = coalesce(family_bp, 0L), family_frac = family_bp / L) %>%
  left_join(sense_core %>% group_by(id) %>%
              summarise(sense_start = min(query_start) - 1L, te_end = max(query_end),
                        retained = sum(cons_end - cons_begin + 1L),
                        divergence = sum(divergence * query_len) / sum(query_len)), by = "id") %>%
  left_join(sense_core %>% group_by(id) %>% slice(which.min(query_start)) %>% ungroup() %>%
              transmute(id, cons_start = cons_begin, five_hit = repeat_name, five_cend = cons_end), by = "id") %>%
  left_join(sense_core %>% group_by(id) %>% slice(which.max(query_end)) %>% ungroup() %>%
              transmute(id, cons_end, cons_left_3p = cons_left, cons_len = cons_end + cons_left, three_hit = repeat_name),
            by = "id") %>%
  left_join(bind_rows(sense_core, hexamer_hits, inverted_core) %>% group_by(id) %>%
              summarise(te_start = min(query_start) - 1L), by = "id") %>%
  mutate(has_main = !is.na(te_end))

## RepeatMasker often leaves the outer copies of the SVA hexamer unannotated or
## calls them CT-rich low complexity, so the element would start partway into
## the array and the rest would be called a 5' transduction. The start is moved
## out over any run of >= 2 hexamer motifs (gaps <= 12 bp) running into it.
## Motifs are every rotation of the canonical CCCTCT and of the alternate units
## of Kirby et al. 2025 (Mobile DNA 16:42, Fig. 2d); their 16/19/25-mer units
## are hexamer copies with short insertions, which the gap tolerance covers.
HEX_UNITS <- c("CCCTCT", "TCTACC", "CCCCCT", "GTCTCC", "CCCCTCT", "CCCCCTCT")
hex_rotations <- unique(unlist(map(HEX_UNITS, function(unit)
  map_chr(seq_len(nchar(unit)), ~ paste0(substring(unit, .x), substr(unit, 1, .x - 1))))))
HEX_MOTIF_RE <- paste(hex_rotations[order(-nchar(hex_rotations))], collapse = "|")
extend_over_hexamer <- calls$family == "SVA" & calls$has_main & calls$te_start > 0
extend_over_hexamer[is.na(extend_over_hexamer)] <- FALSE
extended_start <- calls$te_start
extended_start[extend_over_hexamer] <- map2_int(calls$sense[extend_over_hexamer], calls$te_start[extend_over_hexamer],
                                                function(sense, start) {
  motifs <- stri_locate_all_regex(stri_sub(sense, 1L, start), HEX_MOTIF_RE, omit_no_match = TRUE)[[1]]
  edge <- start; n_motifs <- 0L
  for (k in rev(seq_len(nrow(motifs)))) {
    if (edge - motifs[k, 2] > 12L) break
    edge <- motifs[k, 1] - 1L; n_motifs <- n_motifs + 1L
  }
  if (n_motifs >= 2L) as.integer(edge) else start
})
calls <- calls %>% mutate(hexamer_extended_bp = te_start - extended_start, te_start = extended_start)

## SVA_F1 carries MAST2 exon 1, transcribed from a start inside the exon and
## spliced into the SVA Alu-like region (Kirby et al. 2025 count it as part of
## F1, not a 5' transduction). The exon is GRCh38 chr1:45803532-45804072
## (Ensembl MAST2 exon 1, + strand, the element-sense orientation of F1). In
## the 2026-11 cohort F1 elements start at exon offset 178 and leave the exon at
## 398 or 451, joining the SVA consensus at ~269 or ~392. An SVA with
## f1_min_mast2_bp of MAST2 20-mers in front of its element start is F1; the
## element then starts at the MAST2 start.
MAST2_EXON1 <- paste0(
  "GGGAGGCGGCGCGCGGCCGGAGGTGGCGGGCGGGCTCCCGGTGCGCGCGCGTGAGGGGCGGGGCCGGCCG",
  "AAGGCTGGCTGTAGGCAGGCGGCTGAGCCGGCGGCGGGTGGCCTGCCCAACGTGTGCTGGGTGGGAGAAG",
  "GCGAGGCGTCAGCGATGCTGTCTCTTCCGTGAGGAGCGCAGAGGAGGTCGCGGCGCCGGAGGCCCCAGAA",
  "GGCTCGAAGGCGCCGCGGGCTGGGGTCGGTGGCTTAGGGAGCCCGTCCGGCCATGGTGGCCGCGGGTGGT",
  "GGTTGGCGCGGCTGCGCTGCGGCCCGGGGCAGTGCGGAGCCGGGACAGTCGCGGCGCTGACGCCCGCGGG",
  "CCCCAGCTGCAGATATGAAGCGGAGCCGCTGCCGCGACCGACCGCAGCCGCCGCCGCCCGACCGCCGGGA",
  "GGATGGAGTTCAGCGGGCAGCGGAGCTGTCTCAGTCTTTGCCGCCGCGCCGGCGAGCGCCGCCCGGGAGG",
  "CAGCGGCTGGAGGAGCGGACGGGCCCCGCGGGGCCCGAGGGCAAGGAGCAG")
MAST2_F1_START <- 178L
mast2_kmers <- stri_sub(MAST2_EXON1, seq_len(nchar(MAST2_EXON1) - 19L), length = 20L)
f1_candidates <- calls %>%
  filter(family == "SVA", has_main, te_start >= opt$f1_min_mast2_bp) %>%
  select(id, sense, te_start)
mast2_match <- map2(f1_candidates$sense, f1_candidates$te_start, function(sense, start) {
  front <- stri_sub(sense, 1L, start)
  n_kmers <- nchar(front) - 19L
  at <- if (n_kmers > 0) which(stri_sub(front, seq_len(n_kmers), length = 20L) %in% mast2_kmers) else integer()
  if (!length(at)) return(c(0L, NA_integer_, NA_integer_))
  c(length(unique(unlist(map(at, ~ .x:(.x + 19L))))), min(at), match(stri_sub(front, min(at), length = 20L), mast2_kmers))
})
f1_found <- f1_candidates %>%
  transmute(id, mast2_bp = map_int(mast2_match, 1), mast2_q_start = map_int(mast2_match, 2),
            mast2_exon_start = map_int(mast2_match, 3)) %>%
  filter(mast2_bp >= opt$f1_min_mast2_bp)
calls <- calls %>%
  mutate(subfamily_input = subfamily) %>%
  left_join(f1_found %>% rename(f1_bp = mast2_bp, f1_q_start = mast2_q_start, f1_exon_start = mast2_exon_start), by = "id") %>%
  mutate(sva_f1 = !is.na(f1_bp),
         mast2_bp = coalesce(f1_bp, 0L),
         mast2_q_start = f1_q_start,
         mast2_exon_start = f1_exon_start,
         te_start = if_else(sva_f1, f1_q_start - 1L, te_start),
         subfamily = if_else(sva_f1, "SVA_F1", subfamily)) %>%
  select(-f1_bp, -f1_q_start, -f1_exon_start)

## ---- Classification: truncation, inversion, deletions (Fig. 2 b/c/e/i) ---------------

SPLIT_3P_RE <- "_(3end|orf2)$"
## The SVA consensus starts with ~56 bp of hexamer, so when RepeatMasker gives
## the element's hexamer to a simple-repeat hit the SVA hit starts at
## consensus ~57 and a complete element reads as 5' truncated. Hexamer bases
## directly in front count toward the consensus 5' end, but only while the
## hit starts inside the consensus hexamer (sva_hexamer_cons_bp).
calls <- calls %>%
  mutate(max5 = unname(max_5p_missing[family]), max3 = unname(max_3p_missing[family]),
         hex_before_bp = if_else(family == "SVA" & has_main & !sva_f1 & !id %in% inverted_core$id,
                                 pmax(0L, sense_start - te_start), 0L),
         five_prime_missing_rm = cons_start - 1L,
         five_prime_missing = case_when(
           sva_f1 ~ pmax(0L, mast2_exon_start - MAST2_F1_START),
           (hex_before_bp >= 12L & five_prime_missing_rm <= opt$sva_hexamer_cons_bp) %in% TRUE ~ pmax(0L, five_prime_missing_rm - hex_before_bp),
           TRUE ~ five_prime_missing_rm),
         five_prime_truncated = has_main & (stri_detect_regex(five_hit, SPLIT_3P_RE, case_insensitive = TRUE) |
                                              five_prime_missing > max5),
         retained_cons_frac = pmin(1, retained / cons_len),
         three_prime_missing = cons_left_3p,
         three_prime_truncated = has_main & cons_left_3p > max3 & !stri_detect_regex(three_hit, "_5end$"))

## 5' inversion: opposite-strand family sequence 5' of the sense portion.
## junction > 0: deletion at the inversion point, < 0: duplication (Ostertag & Kazazian 2001)
inversions <- inverted_core %>%
  inner_join(calls %>% filter(has_main) %>% select(id, sense_start), by = "id") %>%
  filter(query_end <= sense_start + opt$inv_slack) %>%
  group_by(id) %>%
  summarise(inversion_len = sum(query_len), inversion_cons_start = min(cons_begin), inversion_cons_end = max(cons_end),
            inner_cend = cons_end[which.max(query_end)])
calls <- calls %>%
  left_join(inversions, by = "id") %>%
  mutate(inversion_len = coalesce(inversion_len, 0L),
         five_prime_inversion = inversion_len >= opt$min_inversion,
         inversion_junction = if_else(five_prime_inversion, cons_start - inner_cend - 1L, NA_integer_)) %>%
  select(-inner_cend)

## Internal deletions / rearrangements between consecutive sense hits of one element
deletions <- sense_core %>%
  group_by(id) %>%
  mutate(next_cons_begin = lead(cons_begin), next_query_start = lead(query_start), next_name = lead(repeat_name)) %>%
  ungroup() %>%
  filter(!is.na(next_cons_begin),
         repeat_name == next_name | stri_detect_regex(next_name, SPLIT_3P_RE) | stri_detect_regex(repeat_name, "_5end$")) %>%
  mutate(cons_gap = next_cons_begin - cons_end - 1L, query_gap = next_query_start - query_end - 1L,
         missing = cons_gap - pmax(query_gap, 0L)) %>%
  group_by(id) %>%
  summarise(internal_deletion_bp = sum(pmax(0L, missing) * (missing >= opt$min_internal_gap)),
            internal_rearrangement = any(cons_gap < -opt$min_internal_gap))
calls <- calls %>%
  left_join(deletions, by = "id") %>%
  mutate(internal_deletion_bp = coalesce(internal_deletion_bp, 0L),
         internal_deletion = internal_deletion_bp > 0,
         internal_rearrangement = coalesce(internal_rearrangement, FALSE))

## snRNA chimeras and other repeats inside the insertion
SNRNA_RE <- "snRNA|scRNA|srpRNA|7SL|U6|U3|7SK"
snrna_hits <- nontail_hits %>%
  filter(stri_detect_regex(repeat_class, SNRNA_RE) | stri_detect_regex(repeat_name, SNRNA_RE))
calls <- calls %>%
  left_join(snrna_hits %>% group_by(id) %>% summarise(snRNA_name = paste(sort(unique(repeat_name)), collapse = ",")),
            by = "id") %>%
  mutate(snRNA_chimera = !is.na(snRNA_name) & has_main) %>%
  left_join(nontail_hits %>%
              filter(!(id %in% snrna_hits$id & repeat_name %in% snrna_hits$repeat_name),
                     !repeat_class %in% c("Simple_repeat", "Low_complexity")) %>%
              group_by(id) %>% summarise(other_repeat_bp = sum(query_len)), by = "id") %>%
  mutate(other_repeat_bp = coalesce(other_repeat_bp, 0L))

## ---- Classification: poly(A) tails ----------------------------------------------------

## Maximal A-rich segments (0-based half-open, absolute) in seqs[from+1 ..]:
## runs of A merged across <= 2 non-A bases; kept if >= min_len and >= frac A.
## Leading sparse runs are trimmed off so they cannot sink a real tail.
find_tails <- function(ids, seqs, from, min_len, frac) {
  window <- stri_sub(seqs, from + 1L)
  runs <- stri_locate_all_regex(window, "A+(?:[^A]{1,2}A+)*", omit_no_match = TRUE)
  n_runs <- map_int(runs, nrow)
  if (!sum(n_runs)) return(tibble(id = character(), s = integer(), e = integer()))
  at <- do.call(rbind, runs)
  run_id <- rep(ids, n_runs); run_from <- rep(from, n_runs)
  run_start <- at[, 1]; run_end <- at[, 2]
  run_seq <- stri_sub(rep(window, n_runs), run_start, run_end)
  long <- nchar(run_seq) >= min_len
  run_id <- run_id[long]; run_from <- run_from[long]; run_start <- run_start[long]; run_end <- run_end[long]
  run_seq <- run_seq[long]
  is_tail <- function(x) nchar(x) >= min_len & stri_count_fixed(x, "A") / nchar(x) >= frac
  passes <- is_tail(run_seq)
  for (iteration in 1:100) {
    todo <- which(!passes & stri_detect_regex(run_seq, "^A+[^A]{1,2}A"))
    if (!length(todo)) break
    trimmed <- stri_replace_first_regex(run_seq[todo], "^A+[^A]{1,2}", "")
    run_start[todo] <- run_start[todo] + (nchar(run_seq[todo]) - nchar(trimmed))
    run_seq[todo] <- trimmed
    passes[todo] <- is_tail(trimmed)
  }
  tibble(id = run_id[passes], s = run_from[passes] + run_start[passes] - 1L, e = run_from[passes] + run_end[passes])
}

## Search from just before the element's 3' end (or the whole sequence when
## there is no sense element); RepeatMasker's own tail hits count too.
calls <- calls %>% mutate(tail_search_from = if_else(has_main, pmax(0L, te_end - 40L), 0L))
tails <- bind_rows(find_tails(calls$id, calls$sense, calls$tail_search_from, opt$min_polya, opt$polya_min_frac),
                   tail_hits %>% transmute(id, s = query_start - 1L, e = query_end)) %>%
  distinct() %>%
  inner_join(calls %>% select(id, L, te_end, tail_search_from), by = "id") %>%
  filter(e > tail_search_from)
terminal_tails <- tails %>%
  filter(e >= L - opt$terminal_slack) %>%
  arrange(id, desc(e - s)) %>%
  distinct(id, .keep_all = TRUE) %>%
  select(id, term_s = s, term_e = e)
first_tails <- tails %>%
  filter(!is.na(te_end), s <= te_end + opt$polya_gap, e > te_end - 10L) %>%
  arrange(id, s) %>%
  distinct(id, .keep_all = TRUE) %>%
  select(id, first_s = s, first_e = e)
tails_after_element <- tails %>%
  filter(!is.na(te_end), s <= te_end + opt$polya_gap, e >= te_end) %>%
  arrange(id, s) %>%
  distinct(id, .keep_all = TRUE) %>%
  select(id, after_s = s)
calls <- calls %>%
  left_join(terminal_tails, by = "id") %>%
  left_join(first_tails, by = "id") %>%
  left_join(tails_after_element, by = "id") %>%
  mutate(terminal_polyA = !is.na(term_s),
         polyA_len = if_else(terminal_polyA, term_e - term_s, 0L),
         internal_polyA = has_main & !is.na(first_s) & terminal_polyA & first_e < term_s)

## ---- Classification: transductions and orphans (Fig. 2 d/f/h) -----------------------

## Fraction of each interval [s, e) (0-based) covered by hits (query_start..query_end, 1-based).
masked_frac <- function(intervals, hits) {
  if (!nrow(intervals)) return(numeric())
  covered <- intervals %>%
    mutate(row = row_number()) %>%
    inner_join(hits %>% select(id, query_start, query_end), by = "id", relationship = "many-to-many") %>%
    group_by(row) %>%
    summarise(covered = sum(pmax(0L, pmin(e, query_end) - pmax(s, query_start - 1L))))
  bp <- numeric(nrow(intervals))
  bp[covered$row] <- covered$covered
  pmin(1, bp / pmax(1, intervals$e - intervals$s))
}

## Simple-repeat / low-complexity sequence at either end of a transduction
## candidate cannot be traced and is usually element (the SVA hexamer called
## CT-rich) or calling noise, so it is trimmed off; interior copies sit
## between unique flank and stay. Intervals are 0-based half-open.
simple_by_id <- split(nontail_hits %>% filter(stri_detect_regex(repeat_class, "Simple_repeat|Low_complexity")) %>%
                        select(id, query_start, query_end), ~ id)
trim_simple <- function(intervals) {
  for (i in seq_len(nrow(intervals))) {
    simple <- simple_by_id[[intervals$id[i]]]
    if (is.null(simple)) next
    s <- intervals$s[i]; e <- intervals$e[i]
    repeat {
      at_start <- simple$query_end[simple$query_start - 1L <= s + 2L & simple$query_end > s & simple$query_start - 1L < e]
      at_end <- simple$query_start[simple$query_end >= e - 2L & simple$query_start - 1L < e & simple$query_end > s]
      if (!length(at_start) && !length(at_end)) break
      if (length(at_start)) s <- max(at_start)
      if (length(at_end)) e <- min(at_end) - 1L
      if (e <= s) break
    }
    intervals$s[i] <- as.integer(s); intervals$e[i] <- as.integer(max(s, e))
  }
  intervals
}

## 3' transduction: unique sequence between the element (or its first, internal
## poly(A)) and the terminal poly(A)
candidates_3p <- calls %>%
  filter(has_main, terminal_polyA) %>%
  transmute(id, s = if_else(internal_polyA, first_e, te_end), e = term_s)
trimmed_3p <- trim_simple(candidates_3p %>% filter(e > s)) %>% rename(trim_s = s, trim_e = e)
candidates_3p <- candidates_3p %>%
  left_join(trimmed_3p, by = "id") %>%
  mutate(trimmed_bp = if_else(is.na(trim_s), 0L, (e - s) - (trim_e - trim_s)),
         s = coalesce(trim_s, s), e = coalesce(trim_e, e)) %>%
  select(-trim_s, -trim_e)
long_3p <- candidates_3p %>% filter(e - s >= opt$min_transduction)
long_3p$masked <- masked_frac(long_3p, nontail_hits)
calls <- calls %>%
  left_join(candidates_3p %>% rename(c3_s = s, c3_e = e, c3_trimmed = trimmed_bp), by = "id") %>%
  left_join(long_3p %>% select(id, tdx3_masked_frac = masked), by = "id") %>%
  mutate(three_prime_transduction = !is.na(tdx3_masked_frac) & tdx3_masked_frac <= opt$max_transduction_masked,
         tdx3_start = if_else(three_prime_transduction, c3_s, NA_integer_),
         tdx3_end = if_else(three_prime_transduction, c3_e, NA_integer_),
         tdx3_len = tdx3_end - tdx3_start,
         tdx3_trimmed_bp = if_else(three_prime_transduction, c3_trimmed, 0L),
         tdx3_type = case_when(three_prime_transduction & internal_polyA ~ "TE-A-unique-A",
                               three_prime_transduction ~ "TE-unique-A",
                               TRUE ~ NA_character_),
         unresolved_3p_extension = has_main & !terminal_polyA & (L - te_end - opt$terminal_slack) >= opt$min_transduction) %>%
  select(-c3_s, -c3_e, -c3_trimmed)

## orphan candidates: no sense element, unique sequence + terminal poly(A)
orphans <- calls %>% filter(!has_main, terminal_polyA) %>% transmute(id, s = 0L, e = term_s)
orphans$unique_bp <- (orphans$e - orphans$s) * (1 - masked_frac(orphans, other_hits))
calls <- calls %>%
  left_join(orphans %>% select(id, orphan_unique_bp = unique_bp), by = "id") %>%
  mutate(orphan_candidate = (orphan_unique_bp >= opt$min_orphan) %in% TRUE,
         tdx3_start = if_else(orphan_candidate, 0L, tdx3_start),
         tdx3_end = if_else(orphan_candidate, term_s, tdx3_end),
         tdx3_len = if_else(orphan_candidate, term_s, tdx3_len)) %>%
  select(-orphan_unique_bp)

## 5' foreign sequence -> 5' transduction (intact element 5' end) or
## exon-trapping candidate (internally starting element, spliced junction)
trimmed_5p <- trim_simple(calls %>% filter(has_main, te_start > 0) %>% transmute(id, s = 0L, e = te_start))
calls <- calls %>%
  left_join(trimmed_5p %>% rename(t5_s = s, t5_e = e), by = "id") %>%
  mutate(tdx5_trimmed_bp = if_else(is.na(t5_s), 0L, te_start - (t5_e - t5_s)))
candidates_5p <- calls %>%
  filter(has_main, !is.na(t5_e), (t5_e - pmax(t5_s, opt$terminal_slack)) >= opt$min_transduction) %>%
  transmute(id, s = pmax(t5_s, opt$terminal_slack), e = t5_e)
candidates_5p$masked <- masked_frac(candidates_5p, nontail_hits)
calls <- calls %>%
  left_join(candidates_5p %>% select(id, tdx5_masked_frac = masked), by = "id") %>%
  mutate(is_candidate_5p = id %in% candidates_5p$id,
         tdx5_t_rich_head = if_else(is_candidate_5p, stri_count_fixed(stri_sub(sense, 1, 20), "T") / 20 >= 0.8, NA),
         foreign_5p = (is_candidate_5p & tdx5_masked_frac <= opt$max_transduction_masked & !tdx5_t_rich_head) %in% TRUE,
         exon_trapping_candidate = foreign_5p & five_prime_truncated & !five_prime_inversion & !sva_f1,
         five_prime_transduction = foreign_5p & !exon_trapping_candidate,
         tdx5_junction_dinuc = if_else(exon_trapping_candidate, stri_sub(sense, pmax(1L, te_start - 1L), te_start), NA_character_),
         called_5p = five_prime_transduction | exon_trapping_candidate,
         tdx5_start = if_else(called_5p, t5_s, NA_integer_),
         tdx5_end = if_else(called_5p, t5_e, NA_integer_),
         tdx5_len = tdx5_end - tdx5_start,
         tdx5_trimmed_bp = if_else(called_5p, tdx5_trimmed_bp, 0L)) %>%
  select(-t5_s, -t5_e, -is_candidate_5p, -foreign_5p, -called_5p)

## ---- Classification: mechanism of 3' truncation and the structure label (Fig. 2 g/i) -----

PAS_RE <- "AATAAA|ATTAAA|AGTAAA|TATAAA|CATAAA|GATAAA|AATATA|AATACA|AATAGA|ACTAAA|AAGAAA|AATGAA"
STRUCTURE_FLAGS <- c("five_prime_truncated", "five_prime_inversion", "three_prime_truncated", "internal_priming",
                     "premature_polyadenylation", "internal_deletion", "internal_rearrangement",
                     "three_prime_transduction", "five_prime_transduction", "exon_trapping_candidate", "snRNA_chimera")
calls <- calls %>%
  mutate(element_end_20 = stri_sub(sense, pmax(1L, te_end - 19L), te_end),
         three_prime_end_A_frac = if_else(three_prime_truncated %in% TRUE,
                                          stri_count_fixed(element_end_20, "A") / nchar(element_end_20), NA_real_),
         pas_anchor = if_else(is.na(after_s), te_end, after_s),
         pas_upstream_of_tail = if_else(three_prime_truncated %in% TRUE,
                                        stri_detect_regex(stri_sub(sense, pmax(1L, pas_anchor - 39L), pmax(0L, pas_anchor - 5L)), PAS_RE),
                                        NA),
         internal_priming = three_prime_truncated & !is.na(three_prime_end_A_frac) & three_prime_end_A_frac >= 0.7,
         premature_polyadenylation = three_prime_truncated & !internal_priming & !is.na(after_s) & pas_upstream_of_tail %in% TRUE,
         full_length = has_main & !five_prime_truncated & !five_prime_inversion & !three_prime_truncated) %>%
  select(-element_end_20, -pas_anchor)
structure_parts <- map(STRUCTURE_FLAGS, ~ if_else(calls[[.x]] %in% TRUE, .x, NA_character_))
structure_parts <- c(structure_parts, list(if_else(calls$has_main & !calls$terminal_polyA, "no_terminal_polyA", NA_character_)))
calls <- calls %>%
  mutate(structure = pmap_chr(structure_parts, function(...) paste(na.omit(c(...)), collapse = "+")),
         structure = case_when(!has_main ~ "no_element", structure == "" ~ "full_length_canonical", TRUE ~ structure),
         tdx3_seq = if_else(!is.na(tdx3_len), stri_sub(sense, tdx3_start + 1L, tdx3_end), NA_character_),
         tdx5_seq = if_else(five_prime_transduction | exon_trapping_candidate, stri_sub(sense, tdx5_start + 1L, tdx5_end), NA_character_)) %>%
  select(-max5, -max3, -tail_search_from, -term_s, -term_e, -first_s, -first_e, -after_s, -sense, -five_cend)

## ---- Per-sequence table -------------------------------------------------------------

FLAGS <- c("full_length", "five_prime_truncated", "three_prime_transduction", "five_prime_inversion",
           "five_prime_transduction", "internal_priming", "exon_trapping_candidate", "three_prime_truncated",
           "premature_polyadenylation", "internal_deletion", "internal_rearrangement", "orphan_candidate",
           "snRNA_chimera", "terminal_polyA", "internal_polyA")
per_seq <- insertions %>%
  select(id, limeaid_id, variant_id, chrom, pos, ALT_md5, uid) %>%
  mutate(seq_md5 = uid) %>%
  inner_join(representatives %>% inner_join(calls, by = "id") %>% select(-id), by = "uid") %>%
  arrange(uid) %>%
  select(-uid) %>%
  mutate(across(all_of(FLAGS), ~ coalesce(.x, FALSE)))

## Genotypes: haplotypes carrying each catalogued allele (one row per allele and
## individual, keyed like the allele table, POS = its `start`). hap1/hap2 say
## which haplotypes carry it; GT is used only where they are missing. A
## variant's AF is its alleles' carried haplotypes over AN = 2 x cohort size.
genotype_raw <- read_tsv(params$genotype_table, col_types = cols(.default = col_character()), progress = FALSE)
if (!all(c("CHROM", "POS", "ALT_md5", "indiv") %in% names(genotype_raw))) {
  stop(params$genotype_table, ": genotype table needs columns CHROM, POS, ALT_md5, indiv")
}
as_haplotype <- function(x) {
  x <- toupper(stri_trim_both(x))
  case_when(x %in% c("TRUE", "T", "1") ~ TRUE, x %in% c("FALSE", "F", "0") ~ FALSE, TRUE ~ NA)
}
gt_alleles <- if ("GT" %in% names(genotype_raw)) stri_split_regex(genotype_raw$GT, "[/|]", simplify = TRUE) else
  matrix(NA_character_, nrow(genotype_raw), 2)
genotypes <- tibble(
  chrom = genotype_raw$CHROM, pos = as.integer(genotype_raw$POS), ALT_md5 = tolower(genotype_raw$ALT_md5),
  indiv = genotype_raw$indiv,
  hap1 = if ("hap1" %in% names(genotype_raw)) as_haplotype(genotype_raw$hap1) else NA,
  hap2 = if ("hap2" %in% names(genotype_raw)) as_haplotype(genotype_raw$hap2) else NA,
  gt1 = gt_alleles[, 1] %in% "1",
  gt2 = (if (ncol(gt_alleles) > 1) gt_alleles[, 2] else NA_character_) %in% "1"
) %>%
  mutate(hap1 = coalesce(hap1, gt1), hap2 = coalesce(hap2, gt2)) %>%
  group_by(chrom, pos, ALT_md5, indiv) %>%      # a repeated row for one allele and individual counts once
  summarise(hap1 = any(hap1), hap2 = any(hap2), .groups = "drop")
n_individuals <- as.integer(params$n_individuals)
if (is.na(n_individuals) || n_individuals <= 0) n_individuals <- n_distinct(genotypes$indiv)
allele_number <- 2L * n_individuals
carriers <- genotypes %>% filter(hap1 | hap2)
allele_counts <- carriers %>%
  group_by(chrom, pos, ALT_md5) %>%
  summarise(allele_count = sum(hap1) + sum(hap2), n_carriers = n(), .groups = "drop")
per_seq <- per_seq %>%
  left_join(allele_counts, by = c("chrom", "pos", "ALT_md5")) %>%
  mutate(allele_count = coalesce(allele_count, 0L), n_carriers = coalesce(n_carriers, 0L),
         af = allele_count / allele_number)
variant_genotypes <- carriers %>%
  inner_join(per_seq %>% distinct(chrom, pos, ALT_md5, variant_id), by = c("chrom", "pos", "ALT_md5"),
             relationship = "many-to-many") %>%
  group_by(variant_id) %>%
  summarise(n_carriers = n_distinct(indiv), allele_count = sum(hap1) + sum(hap2))
n_unmatched <- nrow(anti_join(allele_counts, per_seq, by = c("chrom", "pos", "ALT_md5")))
message(sprintf("genotype table: %d individuals (AN = %d); %d of %d catalogued alleles carried by someone; %d genotyped alleles not in the catalog.",
                n_individuals, allele_number, sum(per_seq$allele_count > 0), nrow(per_seq), n_unmatched))
if (!any(per_seq$allele_count > 0)) stop("No catalogued allele matched the genotype table on CHROM + POS + ALT_md5")
write_tsv(per_seq %>% select(-tdx3_seq, -tdx5_seq), out_file("per_sequence.tsv.gz"), na = "NA")
message("wrote ", out_file("per_sequence.tsv.gz"), " (", nrow(per_seq), " rows)")

## ---- Per-variant rollup -------------------------------------------------------------
## Modal structure, majority flags, medians and concordance, each sequence
## weighted by the haplotypes carrying it (a variant none of whose alleles is
## carried falls back to equal weights).

weighted_mode <- function(x, w) { totals <- tapply(w, x, sum); names(totals)[which.max(totals)] }
weighted_median <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  ord <- order(x)
  as.numeric(x[ord][which(cumsum(w[ord]) >= sum(w) / 2)[1]])
}
MEDIAN_COLS <- c("L", "five_prime_missing", "three_prime_missing", "retained_cons_frac", "polyA_len",
                 "tdx3_len", "tdx5_len", "inversion_len", "divergence")
weighted <- per_seq %>%
  group_by(variant_id) %>%
  mutate(w = as.numeric(allele_count), w = if (sum(w) > 0) w else rep(1, n())) %>%
  ungroup()
per_var <- weighted %>%
  group_by(variant_id) %>%
  summarise(family = family[1], subfamily = weighted_mode(coalesce(subfamily, "NA"), w),
            chrom = chrom[1], pos = pos[1], n_sequences = n(), n_distinct_seqs = n_distinct(seq_md5),
            structure = weighted_mode(structure, w),
            across(all_of(FLAGS), ~ sum(w * .x) / sum(w) >= 0.5),
            across(all_of(MEDIAN_COLS), ~ weighted_median(.x, w)),
            .groups = "drop")
concordance <- weighted %>%
  inner_join(per_var %>% select(variant_id, modal = structure), by = "variant_id") %>%
  group_by(variant_id) %>%
  summarise(structure_concordance = sum(w * (structure == modal)) / sum(w),
            allele_concordance = mean(structure == modal))
per_var <- per_var %>%
  inner_join(concordance, by = "variant_id") %>%
  left_join(variant_genotypes, by = "variant_id") %>%
  mutate(n_carriers = coalesce(n_carriers, 0L), allele_count = coalesce(allele_count, 0L),
         af = allele_count / allele_number) %>%
  arrange(variant_id)
write_tsv(per_var, out_file("per_variant.tsv"), na = "NA")
message("wrote ", out_file("per_variant.tsv"), " (", nrow(per_var), " rows)")
message(sprintf("%d sequences -> %d variants.", nrow(per_seq), nrow(per_var)))

## One transduced segment per variant and kind (the one carried by most
## sequences; ties go to the first seen), named variant|kind|family|chrom|pos
## for source tracing.
transduced <- bind_rows(
  per_seq %>% filter(!is.na(tdx3_seq)) %>% transmute(variant_id, kind = "3p", seq = tdx3_seq, family, chrom, pos),
  per_seq %>% filter(!is.na(tdx5_seq)) %>% transmute(variant_id, kind = "5p", seq = tdx5_seq, family, chrom, pos)
) %>%
  mutate(first_seen = row_number()) %>%
  group_by(variant_id, kind, seq, family, chrom, pos) %>%
  summarise(n = n(), first_seen = min(first_seen), .groups = "drop") %>%
  arrange(desc(n), first_seen) %>%
  distinct(variant_id, kind, .keep_all = TRUE) %>%
  mutate(name = paste(variant_id, kind, family, chrom, pos, sep = "|"))
transductions_fa <- out_file("transductions.fa")
writeLines(paste0(">", transduced$name, "\n", transduced$seq), transductions_fa)
message(sprintf("wrote %s (%d segments)", transductions_fa, nrow(transduced)))

if (!run_sources) {
  message("No --reference / --transduction-paf: skipping transduction source tracing.")
  quit(save = "no", status = 0)
}
if (!nrow(transduced)) {
  message("No transduced segments: nothing to trace.")
  quit(save = "no", status = 0)
}

## ---- Source tracing: genome mapping ------------------------------------------------
## Segments are mapped to the reference; hits within self_window of the
## insertion's own site are flanking-sequence capture, not transduction. Each
## remaining hit is assigned to a same-family, same-strand element upstream
## (3') or downstream (5') of it within source_window. Unmapped/self-only
## segments that recur are clustered all-vs-all: shared non-reference sources.

read_paf <- function(path) {
  lines <- readLines(path)
  if (!length(lines)) return(tibble())
  cols <- stri_split_fixed(lines, "\t", n = 13, simplify = TRUE)
  tibble(qname = cols[, 1], qlen = as.integer(cols[, 2]), qstart = as.integer(cols[, 3]), qend = as.integer(cols[, 4]),
         strand = cols[, 5], tname = cols[, 6], tlen = as.integer(cols[, 7]), tstart = as.integer(cols[, 8]),
         tend = as.integer(cols[, 9]), matches = as.integer(cols[, 10]), alen = as.integer(cols[, 11]),
         mapq = as.integer(cols[, 12]), primary = stri_detect_fixed(cols[, 13], "tp:A:P"))
}
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  is_header <- stri_startswith_fixed(lines, ">")
  record <- cumsum(is_header)
  ids <- stri_extract_first_regex(stri_sub(lines[is_header], 2), "^\\S+")
  seqs <- map_chr(split(lines[!is_header], record[!is_header]), paste, collapse = "")
  tibble(id = ids[as.integer(names(seqs))], seq = toupper(unname(seqs)))
}

threads <- params$threads
transduced_seqs <- read_fasta(transductions_fa)
paf_path <- params$transduction_paf
if (!nzchar(paf_path)) {
  paf_path <- out_file("transductions.paf")
  if (!nzchar(Sys.which("minimap2"))) stop("minimap2 not found on PATH")
  message(sprintf("Mapping %d transduced segments with minimap2 -x %s", nrow(transduced_seqs), params$minimap2_preset))
  status <- system2("minimap2", c("-c", "-x", params$minimap2_preset, "--secondary=yes", "-N", "10", "-p", params$secondary_ratio,
                                  "-t", threads, params$reference, transductions_fa), stdout = paf_path, stderr = FALSE)
  if (status != 0) stop("minimap2 failed (exit ", status, ")")
}
paf_hits <- if (file.exists(paf_path)) read_paf(paf_path) else tibble()
transduction_meta <- transduced_seqs %>%
  mutate(parts = stri_split_fixed(id, "|", simplify = TRUE)) %>%
  transmute(transduction = id, variant_id = parts[, 1], kind = parts[, 2], family = parts[, 3], ins_chrom = parts[, 4],
            ins_pos = suppressWarnings(as.integer(parts[, 5])), td_len = nchar(seq), maps_to_insertion_site = FALSE)
best_hits <- tibble(qname = character(), tname = character(), tstart = integer(), tend = integer(), strand = character(),
                    mapq = integer(), identity = numeric(), qcov = numeric(), n_good_hits = integer())
piece_hits <- tibble()
if (nrow(paf_hits)) {
  paf_hits <- paf_hits %>%
    mutate(identity = matches / alen, qcov = (qend - qstart) / qlen) %>%
    filter(identity >= opt$min_identity) %>%
    inner_join(transduction_meta %>% select(qname = transduction, ins_chrom, ins_pos), by = "qname") %>%
    arrange(qname) %>%
    mutate(self_hit = tname == ins_chrom & !is.na(ins_pos) &
             (abs(tstart - ins_pos) < opt$self_window | abs(tend - ins_pos) < opt$self_window))
  piece_hits <- paf_hits %>%
    filter(qend - qstart >= opt$min_transduction) %>%
    mutate(method = "genome", evalue = NA_real_)
  whole_hits <- paf_hits %>% filter(qcov >= opt$min_qcov)
  transduction_meta <- transduction_meta %>%
    mutate(maps_to_insertion_site = transduction %in% whole_hits$qname[whole_hits$self_hit])
  whole_hits <- whole_hits %>% filter(!self_hit)
  if (nrow(whole_hits)) {
    best_hits <- whole_hits %>%
      arrange(qname, desc(primary), desc(mapq), desc(matches)) %>%
      distinct(qname, .keep_all = TRUE) %>%
      left_join(whole_hits %>% count(qname, name = "n_good_hits"), by = "qname")
  }
}
sources <- transduction_meta %>%
  left_join(best_hits %>% select(transduction = qname, src_chrom = tname, src_start = tstart, src_end = tend,
                                 src_strand = strand, mapq, identity, qcov, n_good_hits), by = "transduction") %>%
  arrange(transduction) %>%
  mutate(mapping = case_when(is.na(src_chrom) & maps_to_insertion_site ~ "self_only",
                             is.na(src_chrom) ~ "unmapped",
                             mapq >= opt$min_mapq ~ "unique",
                             TRUE ~ "multi"))

## ---- Source tracing: candidate source elements ---------------------------------------
## Candidates are the reference RepeatMasker elements (FIXED sources), any
## known_sources BEDs, and every cohort insertion (POLYMORPHIC sources). A
## polymorphic source is absent from the reference, but the flank it transduces
## is ordinary reference sequence beside its insertion point, so the segment
## maps there with no reference element present; the cohort insertion is
## entered as a zero-width element at that point, in its own orientation.

read_elements <- function(path, set) {
  if (grepl("\\.out(\\.gz)?$", path)) {                                   # RepeatMasker .out (whitespace-aligned)
    lines <- readLines(if (grepl("\\.gz$", path)) gzfile(path) else path)
    cols <- stri_split_regex(stri_trim_both(lines[stri_detect_regex(lines, "^\\s*\\d")]), "\\s+", simplify = TRUE)
    elements <- tibble(chrom = cols[, 5], start = as.integer(cols[, 6]) - 1L, end = as.integer(cols[, 7]),
                       strand = if_else(cols[, 9] == "+", "+", "-"), name = cols[, 10], rclass = cols[, 11],
                       div = as.numeric(cols[, 2]) / 100)
  } else {
    table <- read_tsv(path, col_names = FALSE, col_types = cols(.default = col_character()), progress = FALSE)
    if (ncol(table) >= 17 && all(table$X10 %in% c("+", "-"))) {             # UCSC rmsk.txt
      elements <- table %>% transmute(chrom = X6, start = as.integer(X7), end = as.integer(X8), strand = X10, name = X11,
                                      rclass = paste0(X12, "/", X13), div = as.numeric(X3) / 1000)
    } else {                                                                 # BED
      elements <- tibble(chrom = table$X1, start = as.integer(table$X2), end = as.integer(table$X3),
                         name = if (ncol(table) > 3) table$X4 else ".", strand = if (ncol(table) > 5) table$X6 else ".",
                         rclass = if (ncol(table) > 6) table$X7 else "")
    }
  }
  elements %>% mutate(set = set) %>% filter(!is.na(start))
}
merge_intervals <- function(intervals) {
  intervals %>%
    arrange(chrom, start) %>%
    group_by(chrom) %>%
    mutate(block = cumsum(c(TRUE, start[-1] > cummax(end)[-n()]))) %>%
    group_by(chrom, block) %>%
    summarise(start = min(start), end = max(end), .groups = "drop") %>%
    select(-block)
}

variant_strand <- per_seq %>%
  filter(strand %in% c("+", "-")) %>%
  mutate(seen = row_number()) %>%
  group_by(variant_id, strand) %>%
  summarise(n = n(), seen = min(seen), .groups = "drop") %>%
  arrange(variant_id, desc(n), seen) %>%
  distinct(variant_id, .keep_all = TRUE)
cohort_elements <- per_var %>%
  filter(!is.na(chrom), !is.na(pos)) %>%
  transmute(chrom, start = as.integer(pos), end = as.integer(pos), name = variant_id, rclass = as.character(family),
            set = "cohort") %>%
  inner_join(variant_strand %>% select(name = variant_id, strand), by = "name") %>%
  arrange(name)
rmsk <- if (nzchar(params$ref_rmsk)) read_elements(params$ref_rmsk, "reference") else NULL
## Genome hits: old repeats (MIR, L2, ...) map uniquely, but young copies, any
## Alu (~1M conserved copies) and low-complexity sequence can sit identically
## beside the wrong element. The flank library is masked for every repeat: in a
## library that small, even an old repeat can have one meaningless copy.
young_repeats <- if (!is.null(rmsk)) merge_intervals(rmsk %>%
  filter(is.na(div) | div < opt$repeat_max_div | stri_detect_regex(rclass, "Simple_repeat|Low_complexity|(^|/)Alu"))) else NULL
all_repeats <- if (!is.null(rmsk)) merge_intervals(rmsk) else NULL
known_source_paths <- Filter(nzchar, stri_split_fixed(params$known_sources, ",")[[1]])
source_candidates <- bind_rows(c(if (!is.null(rmsk)) list(rmsk), map(known_source_paths, ~ read_elements(.x, basename(.x))),
                                 list(cohort_elements)))
message(sprintf("source candidates: %d reference, %d known_sources, %d cohort (polymorphic) insertion sites",
                sum(source_candidates$set == "reference"), sum(!source_candidates$set %in% c("reference", "cohort")),
                sum(source_candidates$set == "cohort")))
FAMILY_RE <- c(Alu = "^SINE/Alu|(^|/)Alu", L1 = "^LINE/L1|(^|/)L1", SVA = "(^|/)SVA")
source_candidates <- source_candidates %>% filter(set != "reference" | stri_detect_regex(rclass, "Alu|L1|SVA"))
candidate_ends <- source_candidates %>%
  transmute(chrom, pos = end, el_start = start, el_end = end, el_strand = strand, el_name = name, el_class = rclass, set)
candidate_starts <- source_candidates %>%
  transmute(chrom, pos = start, el_start = start, el_end = end, el_strand = strand, el_name = name, el_class = rclass, set)
window <- opt$source_window
tolerance <- opt$overlap_tolerance

## hits: key, variant_id, family, kind, src_chrom, src_start, src_end, src_strand
## -> the nearest qualifying element per key. A 3' transduction's source element
## ends just UPSTREAM of the hit (in the hit's genomic orientation); a 5' one
## starts just DOWNSTREAM.
no_source <- tibble(key = character(), variant_id = character(), family = character(), need = character(),
                    src_start = integer(), el_start = integer(), el_end = integer(), el_strand = character(),
                    el_name = character(), el_class = character(), set = character(), dist = integer())
nearest_source <- function(hits) {
  if (!nrow(source_candidates) || !nrow(hits)) return(no_source)
  hits <- hits %>%
    mutate(upstream = (kind == "3p") == (src_strand == "+"),
           need = if_else(src_strand == "+", "+", "-"),
           lo = if_else(upstream, src_start - window, src_end - tolerance),
           hi = if_else(upstream, src_start + tolerance, src_end + window))
  bind_rows(
    hits %>% filter(upstream) %>%
      inner_join(candidate_ends, join_by(src_chrom == chrom, lo <= pos, hi >= pos)) %>%
      mutate(dist = src_start - el_end),
    hits %>% filter(!upstream) %>%
      inner_join(candidate_starts, join_by(src_chrom == chrom, lo <= pos, hi >= pos)) %>%
      mutate(dist = el_start - src_end)
  ) %>%
    select(key, variant_id, family, need, src_start, el_start, el_end, el_strand, el_name, el_class, set, dist) %>%
    filter(el_strand == need | el_strand == ".", dist >= -tolerance, dist <= window,
           !set %in% c("reference", "cohort") | stri_detect_regex(paste0(el_class, "/", el_name), FAMILY_RE[family]),
           !(set == "cohort" & el_name == variant_id)) %>%
    arrange(key, abs(dist), desc(el_end - el_start)) %>%
    distinct(key, .keep_all = TRUE)
}

assigned <- nearest_source(sources %>% filter(!is.na(src_chrom)) %>%
                             transmute(key = transduction, variant_id, family, kind, src_chrom, src_start, src_end, src_strand))
sources <- sources %>%
  left_join(assigned %>% transmute(transduction = key, source_element = el_name, source_class = el_class, source_set = set,
                                   source_start = el_start, source_end = el_end, source_strand = el_strand,
                                   source_distance = as.integer(dist)), by = "transduction") %>%
  mutate(source_element = as.character(source_element), source_class = as.character(source_class),
         source_set = as.character(source_set), source_start = as.integer(source_start),
         source_end = as.integer(source_end), source_strand = as.character(source_strand),
         source_chrom = if_else(!is.na(source_element), src_chrom, NA_character_),
         source_full_length = case_when(
           family == "L1" & source_set %in% "reference" ~ (source_end - source_start) >= opt$full_length_l1_bp,
           source_set %in% "cohort" ~ per_var$full_length[match(source_element, per_var$variant_id)],
           TRUE ~ NA),
         source_af = if_else(source_set %in% "cohort", per_var$af[match(source_element, per_var$variant_id)], NA_real_),
         source_type = case_when(is.na(source_set) ~ NA_character_, source_set == "reference" ~ "fixed", TRUE ~ "polymorphic"))
if (nzchar(params$genes) && any(!is.na(sources$src_chrom))) {
  genes <- read_elements(params$genes, "gene")
  sources <- sources %>%
    left_join(sources %>% filter(!is.na(src_chrom)) %>%
                inner_join(genes, join_by(src_chrom == chrom, src_end >= start, src_start <= end)) %>%
                group_by(transduction) %>%
                summarise(overlapping_genes = paste(sort(unique(name)), collapse = ",")), by = "transduction")
}

## Group by source: the assigned element, else mapped loci within cluster_window
cluster_window <- opt$cluster_window
sources <- sources %>%
  mutate(source_locus = if_else(!is.na(source_element),
                                sprintf("%s:%d-%d(%s)", source_chrom, source_start, source_end, source_strand), NA_character_))
element_less <- sources %>%
  filter(is.na(source_element), !is.na(src_chrom)) %>%
  arrange(src_chrom, src_start) %>%
  mutate(block = cumsum(c(TRUE, src_chrom[-1] != src_chrom[-n()] | src_start[-1] - cummax(src_end)[-n()] > cluster_window))) %>%
  group_by(block) %>%
  mutate(locus = sprintf("%s:%d(no_element)", src_chrom[1], src_start[1])) %>%
  ungroup()
sources <- sources %>%
  left_join(element_less %>% select(transduction, locus), by = "transduction") %>%
  mutate(source_locus = coalesce(source_locus, locus)) %>%
  select(-locus)

## Unmapped / self-only segments that recur: connected components of an
## all-vs-all alignment are shared non-reference sources.
unplaced <- sources %>% filter(mapping %in% c("unmapped", "self_only"))
if (nrow(unplaced) > 1 && opt$run_ava && nzchar(Sys.which("minimap2"))) {
  unplaced_fa <- file.path(tempdir(), "unmapped_transductions.fa")
  writeLines(paste0(">", unplaced$transduction, "\n", transduced_seqs$seq[match(unplaced$transduction, transduced_seqs$id)]),
             unplaced_fa)
  ava_paf <- file.path(tempdir(), "unmapped_ava.paf")
  status <- system2("minimap2", c("-c", "-x", params$minimap2_preset, "--secondary=yes", "-X", "-t", threads,
                                  unplaced_fa, unplaced_fa), stdout = ava_paf, stderr = FALSE)
  if (status != 0) stop("minimap2 failed (exit ", status, ")")
  ava_hits <- read_paf(ava_paf)
  if (nrow(ava_hits)) {
    ava_hits <- ava_hits %>%
      filter(qname != tname, alen > 0, matches / alen >= opt$min_identity,
             (qend - qstart) / pmin(qlen, tlen) >= opt$min_qcov)
    nodes <- unplaced$transduction
    component <- setNames(seq_along(nodes), nodes)
    from <- match(ava_hits$qname, nodes); to <- match(ava_hits$tname, nodes)
    repeat {   # label propagation to connected components
      label <- pmin(component[from], component[to])
      lowest <- tibble(node = c(from, to), label = c(label, label)) %>% group_by(node) %>% summarise(label = min(label))
      updated <- component
      updated[lowest$node] <- pmin(updated[lowest$node], lowest$label)
      if (identical(updated, component)) break
      component <- updated
    }
    sizes <- table(component)
    shared <- names(sizes)[sizes > 1]
    if (length(shared)) {
      cluster <- match(as.character(component), shared)
      in_cluster <- match(nodes[!is.na(cluster)], sources$transduction)
      sources$source_locus[in_cluster] <- sprintf("nonref_cluster_%04d", cluster[!is.na(cluster)])
    }
  }
}

## Segments with no reference locus (unmapped / self-only, including the
## non-reference clusters): is the sequence carried by some OTHER cohort
## insertion allele? Members of the segment's own non-reference cluster are
## excluded, since they share it by construction.
sources$other_alleles <- NA_character_
unplaced <- sources %>% filter(mapping %in% c("unmapped", "self_only"))
if (nrow(unplaced) && nzchar(Sys.which("minimap2"))) {
  alleles_fa <- file.path(tempdir(), "cohort_alleles.fa")
  writeLines(paste0(">", insertions$variant_id, "|", insertions$id, "\n", insertions$seq), alleles_fa)
  unplaced_fa <- file.path(tempdir(), "unplaced_segments.fa")
  writeLines(paste0(">", unplaced$transduction, "\n", transduced_seqs$seq[match(unplaced$transduction, transduced_seqs$id)]),
             unplaced_fa)
  allele_paf <- file.path(tempdir(), "unplaced_vs_alleles.paf")
  status <- system2("minimap2", c("-c", "-x", params$minimap2_preset, "--secondary=yes", "-N", "50", "-t", threads,
                                  alleles_fa, unplaced_fa), stdout = allele_paf, stderr = FALSE)
  if (status != 0) stop("minimap2 failed (exit ", status, ")")
  allele_hits <- read_paf(allele_paf)
  if (nrow(allele_hits)) {
    cluster_mates <- sources %>%
      filter(stri_startswith_fixed(source_locus, "nonref_cluster") %in% TRUE) %>%
      transmute(mate = paste(source_locus, variant_id))
    found_in <- allele_hits %>%
      filter(alen > 0, matches / alen >= opt$min_identity, (qend - qstart) / qlen >= opt$min_qcov) %>%
      mutate(hit_variant = sub("\\|[^|]*$", "", tname)) %>%
      inner_join(sources %>% select(qname = transduction, variant_id, source_locus), by = "qname") %>%
      filter(hit_variant != variant_id, !paste(source_locus, hit_variant) %in% cluster_mates$mate) %>%
      group_by(qname) %>%
      summarise(found = paste(sort(unique(hit_variant)), collapse = ","))
    sources <- sources %>%
      left_join(found_in, by = c("transduction" = "qname")) %>%
      mutate(other_alleles = coalesce(found, other_alleles)) %>%
      select(-found)
  }
}

## ---- Transduction segments --------------------------------------------------------
## A source element can itself carry a transduction, so one transduced
## segment may hold sequence from several source loci ("concatenated"
## transductions, Kirby et al. 2025, Mobile DNA 16:42). Outward from the
## element, the innermost sub-segment is the oldest and the outermost the
## immediate source's flank. Sub-segments tile the transduction: the
## best-scoring query-disjoint chain of hits; an inherited piece split off the
## inner end of a chain hit when a different locus with its own source element
## aligns there (the immediate source is a reference element that still
## carries the inherited sequence, so the whole segment aligns to it; with
## several such generations the inherited pieces split off as a ladder); and
## stretches >= min_segment with no hit, left unaligned. `inner`/`outer` run
## outward from the element (0-based, half-open).

transduction_segments <- function(piece_hits, transduction_meta, transduced_seqs) {
  lengths_by_q <- transduction_meta %>% select(qname = transduction, qlen = td_len, kind)
  segments <- tibble(qname = character(), inner = integer(), outer = integer(), status = character(), hit_id = integer())
  if (nrow(piece_hits)) {
    piece_hits <- piece_hits %>%
      inner_join(transduction_meta %>% select(qname = transduction, variant_id, family, kind), by = "qname") %>%
      arrange(qname) %>%
      mutate(hit_id = row_number())
    placed <- nearest_source(piece_hits %>% filter(!self_hit) %>%
                               transmute(key = as.character(hit_id), variant_id, family, kind, src_chrom = tname,
                                         src_start = tstart, src_end = tend, src_strand = strand))
    piece_hits <- piece_hits %>%
      left_join(placed %>% transmute(hit_id = as.integer(key), source_element = el_name, source_set = set,
                                     source_start = el_start, source_end = el_end, source_strand = el_strand,
                                     source_distance = as.integer(dist)), by = "hit_id") %>%
      mutate(source_element = as.character(source_element), source_set = as.character(source_set),
             source_start = as.integer(source_start), source_end = as.integer(source_end),
             source_strand = as.character(source_strand),
             status = case_when(self_hit ~ "insertion_site",
                                mapq >= opt$min_mapq | !is.na(source_element) ~ "traced",
                                TRUE ~ "multi_mapping"),
             inner = if_else(kind == "3p", qstart, qlen - qend),
             outer = if_else(kind == "3p", qend, qlen - qstart))
    ## a low-MAPQ hit placed only by a nearby element is a guess when an equally
    ## good hit covers the same stretch elsewhere: segmental duplications put
    ## identical flanks beside several elements
    low_mapq <- piece_hits %>% filter(status == "traced", mapq < opt$min_mapq)
    if (nrow(low_mapq)) {
      tied <- low_mapq %>%
        select(qname, hit_id, inner, outer, matches, tname, tstart) %>%
        inner_join(piece_hits %>% select(qname, rival = hit_id, rival_inner = inner, rival_outer = outer,
                                         rival_matches = matches, rival_chrom = tname, rival_start = tstart),
                   by = "qname", relationship = "many-to-many") %>%
        filter(rival != hit_id,
               pmin(outer, rival_outer) - pmax(inner, rival_inner) >= 0.8 * pmax(outer - inner, rival_outer - rival_inner),
               rival_matches >= 0.95 * matches,
               rival_chrom != tname | abs(rival_start - tstart) > cluster_window)
      piece_hits <- piece_hits %>% mutate(status = if_else(hit_id %in% tied$hit_id, "ambiguous_source", status))
    }
    ## genome-wide, minimap2 may not even report the rival copy; a low-MAPQ
    ## genome hit is kept only as a ladder rung (inside a longer hit from the
    ## same inner end), and anything else is left to the flank search, which
    ## sees every candidate copy
    low_mapq <- piece_hits %>% filter(status == "traced", mapq < opt$min_mapq, method %in% "genome")
    if (nrow(low_mapq)) {
      rungs <- low_mapq %>%
        select(qname, hit_id, inner, outer) %>%
        inner_join(piece_hits %>% select(qname, other = hit_id, other_inner = inner, other_outer = outer),
                   by = "qname", relationship = "many-to-many") %>%
        filter(other != hit_id, abs(other_inner - inner) <= pmin(tolerance, 0.2 * (outer - inner)),
               other_outer >= outer + opt$min_segment)
      piece_hits <- piece_hits %>%
        mutate(status = if_else(hit_id %in% setdiff(low_mapq$hit_id, rungs$hit_id), "ambiguous_source", status))
    }
    ## mapping quality says nothing about a hit that is all repeat (an Alu
    ## fragment, SVA VNTR): it must hold min_segment bp of unique sequence
    repeat_overlap <- if (is.null(young_repeats)) tibble(hit_id = integer(), repeat_bp = numeric()) else
      piece_hits %>%
        filter(status == "traced") %>%
        select(hit_id, tname, tstart, tend) %>%
        inner_join(young_repeats, join_by(tname == chrom, tstart < end, tend > start)) %>%
        group_by(hit_id) %>%
        summarise(repeat_bp = as.numeric(sum(pmin(end, tend) - pmax(start, tstart))))
    piece_hits <- piece_hits %>%
      left_join(repeat_overlap, by = "hit_id") %>%
      mutate(unique_bp = if_else(status == "traced", (tend - tstart) - coalesce(repeat_bp, 0), NA_real_),
             status = if_else(status == "traced" & unique_bp < opt$min_segment, "repeat_only", status),
             chain_weight = if_else(status %in% c("multi_mapping", "ambiguous_source", "repeat_only"), 0.1, 1) * matches)
    ## best-scoring chain of query-disjoint hits per transduction (overlap
    ## allowed up to min(tolerance, half the shorter hit))
    chain <- piece_hits %>%
      select(-repeat_bp) %>%
      group_by(qname) %>%
      group_modify(function(hits, key) {
        hits <- hits %>% arrange(outer)
        score <- numeric(nrow(hits)); previous <- integer(nrow(hits))
        for (i in seq_len(nrow(hits))) {
          earlier <- seq_len(i - 1L)
          fits <- earlier[hits$outer[earlier] - hits$inner[i] <=
                            pmin(tolerance, (hits$outer[i] - hits$inner[i]) / 2, (hits$outer[earlier] - hits$inner[earlier]) / 2)]
          best <- if (length(fits)) fits[which.max(score[fits])] else 0L
          score[i] <- hits$chain_weight[i] + if (best) score[best] else 0
          previous[i] <- best
        }
        i <- which.max(score); keep <- integer()
        while (i > 0L) { keep <- c(i, keep); i <- previous[i] }
        hits[keep, ]
      }) %>%
      ungroup() %>%
      mutate(inherited_nested = FALSE)
    ## each older source's copy ends further in than the next one's, so the
    ## hits form a ladder and each rung's outer end is a generation boundary
    nested <- chain %>%
      select(qname, container = hit_id, container_inner = inner, container_outer = outer,
             container_chrom = tname, container_start = tstart) %>%
      inner_join(piece_hits %>% filter(status == "traced", !is.na(source_element), !hit_id %in% chain$hit_id),
                 by = "qname", relationship = "many-to-many") %>%
      filter(abs(inner - container_inner) <= pmin(tolerance, 0.2 * (outer - inner)),
             outer <= container_outer - opt$min_segment,
             tname != container_chrom | abs(tstart - container_start) > cluster_window) %>%
      arrange(container, outer, desc(matches)) %>%
      group_by(container) %>%
      filter(c(TRUE, diff(outer) >= opt$min_segment)) %>%
      mutate(inner = lag(outer, default = container_inner[1])) %>%
      ungroup()
    if (nrow(nested)) {
      chain <- chain %>%
        left_join(nested %>% group_by(container) %>% summarise(top = max(outer)), by = c("hit_id" = "container")) %>%
        mutate(inner = coalesce(top, inner)) %>%
        select(-top)
    }
    segments <- bind_rows(chain, nested %>% mutate(inherited_nested = TRUE) %>%
                            select(-container, -container_inner, -container_outer, -container_chrom, -container_start)) %>%
      arrange(qname, inner) %>%
      group_by(qname) %>%
      mutate(inner = pmax(inner, lag(cummax(outer), default = 0L))) %>%
      ungroup() %>%
      filter(outer - inner >= opt$min_segment)
  }
  gaps <- lengths_by_q %>%
    select(qname, qlen) %>%
    left_join(segments %>% select(qname, inner, outer), by = "qname") %>%
    arrange(qname, inner) %>%
    group_by(qname) %>%
    reframe(gap_inner = c(0L, outer[!is.na(outer)]), gap_outer = c(inner[!is.na(inner)], qlen[1])) %>%
    transmute(qname, inner = gap_inner, outer = gap_outer, status = "unaligned") %>%
    filter(outer - inner >= opt$min_segment)
  segments <- bind_rows(segments, gaps)
  for (column in c("tname", "tstart", "tend", "strand", "mapq", "identity", "matches", "alen", "qstart", "qend",
                   "source_element", "source_set", "source_start", "source_end", "source_strand", "source_distance",
                   "inherited_nested", "method", "evalue")) {
    if (!column %in% names(segments)) segments[[column]] <- NA
  }
  segments <- segments %>%
    select(-any_of(c("kind", "qlen", "variant_id", "family", "ins_chrom", "ins_pos"))) %>%
    inner_join(lengths_by_q %>% inner_join(transduction_meta %>% select(qname = transduction, variant_id, family, ins_chrom, ins_pos),
                                           by = "qname"), by = "qname") %>%
    arrange(qname, inner) %>%
    group_by(qname) %>%
    mutate(segment = row_number()) %>%
    ungroup() %>%
    mutate(seg_start = if_else(kind == "3p", inner + 1L, qlen - outer + 1L),
           seg_end = if_else(kind == "3p", outer, qlen - inner),
           seg_len = seg_end - seg_start + 1L,
           seg_seq = stri_sub(transduced_seqs$seq[match(qname, transduced_seqs$id)], seg_start, seg_end),
           status = if_else(status == "unaligned" & stri_count_fixed(seg_seq, "A") / seg_len >= opt$polya_min_frac, "polyA", status),
           source_locus = case_when(status != "traced" ~ NA_character_,
                                    !is.na(source_element) ~ sprintf("%s:%d-%d(%s)", tname, source_start, source_end, source_strand),
                                    TRUE ~ sprintf("%s:%d(no_element)", tname, tstart)),
           source_type = case_when(is.na(source_set) ~ NA_character_, source_set == "reference" ~ "fixed", TRUE ~ "polymorphic"),
           inherited_nested = coalesce(inherited_nested, FALSE))
  segments$n_other_hits <- 0L
  if (nrow(piece_hits)) {
    overlapping <- segments %>%
      select(qname, segment, seg_inner = inner, seg_outer = outer, hit_id) %>%
      inner_join(piece_hits %>% select(qname, other = hit_id, inner, outer), by = "qname", relationship = "many-to-many") %>%
      filter(is.na(hit_id) | other != hit_id,
             pmin(seg_outer, outer) - pmax(seg_inner, inner) >= 0.5 * (seg_outer - seg_inner)) %>%
      count(qname, segment, name = "n_other")
    segments <- segments %>%
      left_join(overlapping, by = c("qname", "segment")) %>%
      mutate(n_other_hits = coalesce(n_other, 0L)) %>%
      select(-n_other)
  }
  ## 3': a parent's own poly(A) sits between its inherited sequence and its flank
  segments %>%
    mutate(outer_polyA = if_else(kind == "3p",
                                 map2_lgl(qname, seg_end, ~ any(abs(transduction_tails$s[transduction_tails$id == .x] - .y) <= opt$polya_gap)),
                                 NA))
}

transduction_tails <- find_tails(transduced_seqs$id, transduced_seqs$seq, rep(0L, nrow(transduced_seqs)),
                                 opt$min_polya, opt$polya_min_frac)
segments <- transduction_segments(piece_hits, transduction_meta, transduced_seqs)

## ---- Source tracing: second pass against source flanks --------------------------------
## Stretches too short or too repetitive to place genome-wide are searched
## against the flanks of candidate source elements only, a library small enough
## (~10-50 Mb) for a 20-30 bp match to be significant. Every RefRmsk repeat in
## the library and every homopolymer run in the stretches is hard-masked first
## (and BLAST's DUST filter is on): in a small library a repeat fragment can
## have a single, meaningless copy, and poly(A) linkers match A-rich flanks. A
## hit is kept only without a near-tie at another locus, and must pass the same
## source-element rules as a genome hit.

mask_n <- function(seqs, from, to) {
  for (i in which(lengths(from) > 0)) {
    stri_sub_all(seqs[i], from[[i]], to[[i]]) <- strrep("N", to[[i]] - from[[i]] + 1L)
  }
  seqs
}

## Each run of consecutive unplaced segments is one stretch, so a piece is never
## cut short by a neighbouring junk hit.
UNPLACED_STATUS <- c("unaligned", "multi_mapping", "ambiguous_source", "repeat_only", "polyA")
stretches <- segments %>%
  arrange(qname, seg_start) %>%
  group_by(qname) %>%
  mutate(run = cumsum(!status %in% UNPLACED_STATUS | c(TRUE, seg_start[-1] > seg_end[-n()] + 1L))) %>%
  filter(status %in% UNPLACED_STATUS) %>%
  group_by(qname, run) %>%
  summarise(seg_start = min(seg_start), seg_end = max(seg_end), all_polyA = all(status == "polyA"),
            qlen = qlen[1], variant_id = variant_id[1], family = family[1], kind = kind[1],
            ins_chrom = ins_chrom[1], ins_pos = ins_pos[1], .groups = "drop") %>%
  filter(!all_polyA) %>%
  mutate(seg_seq = stri_sub(transduced_seqs$seq[match(qname, transduced_seqs$id)], seg_start, seg_end)) %>%
  filter(seg_end - seg_start + 1L >= opt$min_segment)
flank_hits <- tibble()
search_families <- intersect(intersect(stri_split_fixed(params$flank_families, ",")[[1]], names(FAMILY_RE)),
                             unique(stretches$family))
if (opt$flank_search && nrow(stretches) && length(search_families)) {
  class_name <- paste0(source_candidates$rclass, "/", source_candidates$name)
  in_family <- Reduce(`|`, map(FAMILY_RE[search_families], ~ stri_detect_regex(class_name, .x)))
  short_l1 <- source_candidates$set == "reference" & stri_detect_regex(class_name, FAMILY_RE[["L1"]]) &
    (source_candidates$end - source_candidates$start) < opt$full_length_l1_bp
  flank_sources <- source_candidates[(in_family | !source_candidates$set %in% c("reference", "cohort")) & !short_l1, ]

  work_dir <- file.path(tempdir(), "flank_search"); dir.create(work_dir, showWarnings = FALSE)
  fasta <- file.path(work_dir, "ref.fa")
  if (!file.exists(fasta)) file.symlink(normalizePath(reference_fasta), fasta)
  if (file.exists(paste0(reference_fasta, ".gzi")) && !file.exists(paste0(fasta, ".gzi"))) {
    file.symlink(normalizePath(paste0(reference_fasta, ".gzi")), paste0(fasta, ".gzi"))
  }
  fai <- if (nzchar(params$reference_fasta_index)) params$reference_fasta_index else
    if (file.exists(paste0(reference_fasta, ".fai"))) paste0(reference_fasta, ".fai") else ""
  if (!nzchar(fai)) {
    if (system2("samtools", c("faidx", fasta)) != 0) stop("samtools faidx failed")
  } else if (!file.exists(paste0(fasta, ".fai"))) {
    file.symlink(normalizePath(fai), paste0(fasta, ".fai"))
  }
  contig_lengths <- read_tsv(paste0(fasta, ".fai"), col_names = FALSE,
                             col_types = cols_only(X1 = col_character(), X2 = col_integer()), progress = FALSE) %>%
    rename(chrom = X1, len = X2)

  flank_regions <- bind_rows(flank_sources %>% transmute(chrom, s = start - window, e = start),
                             flank_sources %>% transmute(chrom, s = end, e = end + window)) %>%
    inner_join(contig_lengths, by = "chrom") %>%
    mutate(s = pmax(0, s), e = pmin(len, e)) %>%
    filter(e > s) %>%
    arrange(chrom, s) %>%
    group_by(chrom) %>%
    mutate(block = cumsum(c(TRUE, s[-1] > cummax(e)[-n()]))) %>%
    group_by(chrom, block) %>%
    summarise(s = min(s), e = max(e), .groups = "drop") %>%
    mutate(rid = paste0("r", row_number()), ri = row_number())
  writeLines(sprintf("%s:%.0f-%.0f", flank_regions$chrom, flank_regions$s + 1, flank_regions$e), file.path(work_dir, "regions.txt"))
  if (system2("samtools", c("faidx", fasta, "-r", file.path(work_dir, "regions.txt")),
              stdout = file.path(work_dir, "flanks_raw.fa")) != 0) stop("samtools faidx failed")
  flank_seqs <- read_fasta(file.path(work_dir, "flanks_raw.fa"))$seq
  if (!is.null(all_repeats)) {
    masked <- flank_regions %>%
      select(ri, chrom, s, e) %>%
      inner_join(all_repeats, join_by(chrom, s < end, e > start)) %>%
      transmute(ri, from = as.integer(pmax(start, s) - s + 1), to = as.integer(pmin(end, e) - s)) %>%
      filter(to >= from) %>%
      arrange(ri, from) %>%
      group_by(ri) %>%
      mutate(block = cumsum(from > lag(cummax(to), default = -1L))) %>%
      group_by(ri, block) %>%
      summarise(from = min(from), to = max(to), .groups = "drop")
    by_region <- split(masked, factor(masked$ri, levels = seq_len(nrow(flank_regions))))
    flank_seqs <- mask_n(flank_seqs, map(by_region, "from"), map(by_region, "to"))
  }
  message(sprintf("source-flank search: %.0f%% of the flank library repeat-masked",
                  100 * sum(stri_count_fixed(flank_seqs, "N")) / sum(nchar(flank_seqs))))
  writeLines(paste0(">", flank_regions$rid, "\n", flank_seqs), file.path(work_dir, "flanks.fa"))
  if (system2("makeblastdb", c("-in", file.path(work_dir, "flanks.fa"), "-dbtype", "nucl", "-out", file.path(work_dir, "flanks")),
              stdout = FALSE) != 0) stop("makeblastdb failed")

  stretches <- stretches %>% mutate(qid = paste0("q", row_number()))
  homopolymers <- stri_locate_all_regex(stretches$seg_seq, sprintf("A{%1$d,}|C{%1$d,}|G{%1$d,}|T{%1$d,}", opt$min_polya),
                                        omit_no_match = TRUE)
  stretches$query_seq <- mask_n(stretches$seg_seq, map(homopolymers, ~ .x[, 1]), map(homopolymers, ~ .x[, 2]))
  writeLines(paste0(">", stretches$qid, "\n", stretches$query_seq), file.path(work_dir, "stretches.fa"))
  blast_out <- file.path(work_dir, "stretches.blast.tsv")
  if (system2("blastn", c("-task", "blastn-short", "-query", file.path(work_dir, "stretches.fa"), "-db", file.path(work_dir, "flanks"),
                          "-evalue", params$flank_max_evalue, "-dust", "yes", "-max_target_seqs", "50", "-num_threads", threads,
                          "-outfmt", shQuote("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore nident"),
                          "-out", blast_out)) != 0) stop("blastn failed")
  blast_hits <- if (file.size(blast_out) > 0) {
    read_tsv(blast_out, col_names = c("qid", "rid", "pident", "length", "qs", "qe", "ss", "se", "evalue", "bitscore", "nident"),
             col_types = "ccdiiiiiddi", progress = FALSE)
  } else tibble()
  message(sprintf("source-flank search: %d untraced stretches against %d flank regions (%.1f Mb) of %d candidate sources",
                  nrow(stretches), nrow(flank_regions), sum(flank_regions$e - flank_regions$s) / 1e6, nrow(flank_sources)))
  if (nrow(blast_hits)) {
    ## A stretch can hold several pieces, so every query-disjoint hit is kept,
    ## best first, unless a near-tie covering the same piece lands anywhere else
    ## (another element, or a second spot in the same flank: a tandem repeat).
    flank_hits <- blast_hits %>%
      filter(pident / 100 >= opt$min_identity, length >= opt$min_segment) %>%
      inner_join(flank_regions %>% select(rid, chrom, region_start = s), by = "rid") %>%
      mutate(strand = if_else(ss <= se, "+", "-"),
             tstart = as.integer(region_start + pmin(ss, se) - 1),
             tend = as.integer(region_start + pmax(ss, se))) %>%
      arrange(rid) %>%
      arrange(qid, desc(bitscore)) %>%
      group_by(qid) %>%
      group_modify(function(hits, key) {
        keep <- integer()
        for (i in seq_len(nrow(hits))) {
          overlap <- pmin(hits$qe, hits$qe[i]) - pmax(hits$qs, hits$qs[i]) + 1
          rival <- overlap >= 0.5 * (hits$qe[i] - hits$qs[i] + 1) & hits$bitscore >= hits$bitscore[i] - 2 &
            (hits$chrom != hits$chrom[i] | abs(hits$tstart - hits$tstart[i]) > hits$qe[i] - hits$qs[i] + 1)
          taken <- length(keep) && any(pmin(hits$qe[keep], hits$qe[i]) - pmax(hits$qs[keep], hits$qs[i]) + 1 > tolerance)
          if (!any(rival) && !taken) keep <- c(keep, i)
        }
        hits[keep, ]
      }) %>%
      ungroup() %>%
      inner_join(stretches %>% select(qid, qname, qlen, seg_start, variant_id, family, kind, ins_chrom, ins_pos), by = "qid") %>%
      transmute(qname, qlen, qstart = seg_start - 1L + qs - 1L, qend = seg_start - 1L + qe, strand, tname = chrom,
                tstart, tend, matches = nident, alen = length, mapq = 0L, primary = TRUE, identity = pident / 100,
                ins_chrom, ins_pos, method = "source_flank", evalue, variant_id, family, kind) %>%
      mutate(self_hit = tname == ins_chrom & !is.na(ins_pos) &
               (abs(tstart - ins_pos) < opt$self_window | abs(tend - ins_pos) < opt$self_window),
             key = as.character(row_number()))
    placed <- nearest_source(flank_hits %>% filter(!self_hit) %>%
                               transmute(key, variant_id, family, kind, src_chrom = tname, src_start = tstart,
                                         src_end = tend, src_strand = strand))
    flank_hits <- flank_hits %>% filter(self_hit | key %in% placed$key)
    message(sprintf("source-flank search: %d stretches placed (%d beside a source element, %d at the insertion's own site)",
                    nrow(flank_hits), sum(!flank_hits$self_hit), sum(flank_hits$self_hit)))
    flank_hits <- flank_hits %>% select(-key, -variant_id, -family, -kind)
  }
}
if (nrow(flank_hits)) segments <- transduction_segments(bind_rows(piece_hits, flank_hits), transduction_meta, transduced_seqs)

## ---- Transduction tables --------------------------------------------------------------

## Distinct loci among traced sub-segments: a new locus wherever the
## chromosome changes or two placements are cluster_window apart with
## different source loci.
source_loci <- segments %>%
  filter(status == "traced") %>%
  arrange(qname, tname, tstart) %>%
  group_by(qname) %>%
  summarise(n_source_loci = 1L + sum(tname != lag(tname) |
                                       (tstart - lag(tstart) > cluster_window & source_locus != lag(source_locus)), na.rm = TRUE))
UNTRACED_STATUS <- c("unaligned", "multi_mapping", "ambiguous_source", "repeat_only")
transductions <- segments %>%
  arrange(qname, segment) %>%
  group_by(qname) %>%
  summarise(n_segments = n(),
            segments = paste(sprintf("%d:%d-%d", segment, seg_start, seg_end), collapse = ";"),
            segment_sources = paste(if_else(status == "traced", source_locus, status), collapse = ";"),
            traced_bp = sum(seg_len[status == "traced"]),
            untraced_bp = sum(seg_len[status %in% UNTRACED_STATUS]),
            n_flank_placed = sum(status == "traced" & method %in% "source_flank"),
            has_insertion_site = any(status == "insertion_site"),
            has_inherited_nested = any(inherited_nested),
            outer_polyA_breaks = sum(outer_polyA[segment < n()] %in% TRUE)) %>%
  left_join(source_loci, by = "qname") %>%
  mutate(n_source_loci = coalesce(n_source_loci, 0L),
         flag = case_when(n_source_loci >= 2L ~ "concatenated",
                          n_source_loci == 1L & untraced_bp == 0L ~ "single_source",
                          n_source_loci == 1L ~ "single_source_partial",
                          has_insertion_site ~ "insertion_site_only",
                          TRUE ~ "untraced")) %>%
  inner_join(transduction_meta %>% select(qname = transduction, variant_id, kind, family, ins_chrom, ins_pos, td_len), by = "qname") %>%
  left_join(transduction_tails %>% count(id, name = "tail_runs"), by = c("qname" = "id")) %>%
  mutate(sequence = transduced_seqs$seq[match(qname, transduced_seqs$id)],
         internal_polyA_runs = if_else(kind == "3p", coalesce(tail_runs, 0L), NA_integer_)) %>%
  arrange(qname)
write_tsv(transductions %>%
            select(variant_id, kind, family, ins_chrom, ins_pos, td_len, flag, n_segments, n_source_loci, segments,
                   segment_sources, traced_bp, untraced_bp, n_flank_placed, has_inherited_nested, outer_polyA_breaks,
                   internal_polyA_runs, sequence),
          out_file("transductions.tsv"), na = "NA")
message("wrote ", out_file("transductions.tsv"), " (", nrow(transductions), " rows)")
write_tsv(segments %>%
            arrange(variant_id, kind, segment) %>%
            transmute(variant_id, kind, family, ins_chrom, ins_pos, segment, status, inherited_nested, seg_start, seg_end, seg_len,
                      chrom = tname, t_start = tstart, t_end = tend, strand, mapq, identity, matches, aligned_bp = alen,
                      hit_q_start = qstart + 1L, hit_q_end = qend, method, evalue, n_other_hits, source_element, source_type,
                      source_locus, source_distance, outer_polyA, seg_seq),
          out_file("transduction_segments.tsv"), na = "NA")
message("wrote ", out_file("transduction_segments.tsv"), " (", nrow(segments), " rows)")
message("transductions: ", paste(sprintf("%s %d", names(table(transductions$flag)), table(transductions$flag)), collapse = ", "))

sources <- sources %>%
  mutate(outcome = case_when(stri_startswith_fixed(source_locus, "nonref_cluster") %in% TRUE ~ "shared_unmapped",
                             source_type %in% "fixed" ~ "assigned_fixed",
                             source_type %in% "polymorphic" ~ "assigned_polymorphic",
                             TRUE ~ mapping))
write_tsv(sources, out_file("transduction_sources.tsv"), na = "NA")
message("wrote ", out_file("transduction_sources.tsv"), " (", nrow(sources), " rows)")

## A non-reference cluster whose offspring all sit at one site is more likely a
## single polymorphic site split across variant ids than a shared source.
source_summary <- sources %>%
  filter(!is.na(source_locus)) %>%
  group_by(source_locus) %>%
  summarise(n_offspring = n_distinct(variant_id), kinds = paste(sort(unique(kind)), collapse = ","),
            families = paste(sort(unique(family)), collapse = ","), source_element = source_element[1],
            source_type = source_type[1], source_af = source_af[1], source_full_length = source_full_length[1],
            median_td_len = as.numeric(median(td_len)), offspring_chroms = paste(sort(unique(ins_chrom)), collapse = ","),
            offspring_span_bp = if (n_distinct(ins_chrom) == 1) as.numeric(diff(range(ins_pos))) else NA_real_,
            other_alleles = paste(sort(unique(unlist(strsplit(other_alleles[!is.na(other_alleles)], ",", fixed = TRUE)))),
                                  collapse = ",")) %>%
  mutate(other_alleles = na_if(other_alleles, ""),
         same_site = if_else(stri_startswith_fixed(source_locus, "nonref_cluster"),
                             !is.na(offspring_span_bp) & offspring_span_bp <= cluster_window, NA)) %>%
  arrange(desc(n_offspring))
write_tsv(source_summary, out_file("source_summary.tsv"), na = "NA")
message("wrote ", out_file("source_summary.tsv"), " (", nrow(source_summary), " rows)")
nonref_clusters <- source_summary %>% filter(stri_startswith_fixed(source_locus, "nonref_cluster"))
message(sprintf("sources: %d fixed, %d polymorphic, %d mapped loci without an element; %d non-reference clusters (%d at a single site, %d found in other cohort alleles)",
                sum(source_summary$source_type %in% "fixed"), sum(source_summary$source_type %in% "polymorphic"),
                sum(is.na(source_summary$source_type) & !stri_startswith_fixed(source_summary$source_locus, "nonref_cluster")),
                nrow(nonref_clusters), sum(nonref_clusters$same_site %in% TRUE), sum(!is.na(nonref_clusters$other_alleles))))
