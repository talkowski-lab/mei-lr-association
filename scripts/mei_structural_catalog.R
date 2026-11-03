#!/usr/bin/env Rscript
## Catalogs the structural forms of non-reference Alu/L1/SVA insertions
## (Hancks & Kazazian 2016, Mobile DNA 7:9, Fig. 2 panels b-i) from the
## RepeatMasker hits L1ME-AID stores per sequence, and optionally traces the
## genomic sources of transduced sequence. Writes tables only; figures and
## tests live in notebooks/mei_structural_catalog_analysis.Rmd.
##
## Outputs (<outdir>/<prefix>.structure_*):
##   per_sequence.tsv.gz      one row per input sequence, every call
##   per_variant.tsv          rollup per variant: modal structure, majority flags, medians, concordance
##   transductions.fa         one transduced segment per variant and kind (variant|kind|family|chrom|pos)
##   transduction_sources.tsv source tracing per segment     (only with --reference / --transduction-paf)
##   source_summary.tsv       offspring per source locus     (only with --reference / --transduction-paf)

suppressPackageStartupMessages({
  library(argparser)
  library(data.table)
  library(stringi)
})

## ---- Input formats -----------------------------------------------------------
##
## The catalog starts from one row per insertion sequence plus the
## RepeatMasker hits that annotate it, read from L1ME-AID output TSVs
## (github.com/Markloftus/L1ME-AID) or the INTACT_MEI
## `*.limeaid.v134.intact.tsv` (same table, FLAGS column overwritten with
## INTACT / INTACT_3end), normalised to an `ins` table and a long `hits` table.
## Per-hit query and consensus coordinates exist only in `Element_Hits` (v1.4)
## / `TE_Hits` (older), stored as a Python-list repr of raw RepeatMasker lines
## -- everything downstream is built on those. Columns are looked up by header
## name because their names/order drift between L1ME-AID versions.
##
## Each sequence is joined to the clustered allele table on chrom + POS + md5
## of the sequence, which gives its variant and family (MEI_type), and then to
## the genotype table on that allele's CHROM + POS + ALT_md5.
##
## RepeatMasker coordinates are 1-based inclusive. On 'C'-strand hits the
## consensus columns are "(left) end begin", so they are re-ordered here and
## cbeg <= cend always holds afterwards.

## ---- Option coercion (argparser hands every option over as a string) -------

p_num <- function(x) as.numeric(x)
p_int <- function(x) as.integer(x)
p_bool <- function(x) isTRUE(as.logical(x)) || identical(tolower(as.character(x)), "true")
p_str <- function(x) if (is.null(x) || is.na(x) || identical(x, "")) NULL else as.character(x)

## "Alu=20,L1=50" -> named list
p_map <- function(x) {
  x <- p_str(x)
  if (is.null(x)) return(list())
  kv <- stri_split_fixed(stri_split_fixed(x, ",")[[1]], "=", n = 2)
  setNames(lapply(kv, function(z) stri_trim_both(z[2])), vapply(kv, function(z) stri_trim_both(z[1]), ""))
}

which_tool <- function(...) {
  for (n in c(...)) {
    if (!is.null(n) && nzchar(n) && nzchar(Sys.which(n))) return(unname(Sys.which(n)))
  }
  NULL
}

run_tool <- function(cmd, args, stdout = "", stderr = FALSE) {
  status <- system2(cmd, args, stdout = stdout, stderr = stderr)
  if (!identical(as.integer(status), 0L)) stop(sprintf("%s failed (exit %s)", cmd, status))
  invisible(status)
}

## fread that also reads .gz without needing R.utils
fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...) else fread(path, ...)
}

## ---- Sequence helpers -------------------------------------------------------

revcomp <- function(x) stri_reverse(chartr("ACGTNacgtn", "TGCANtgcan", x))

read_fasta <- function(path) {
  l <- readLines(path, warn = FALSE)
  l <- l[nzchar(l)]
  h <- stri_startswith_fixed(l, ">")
  grp <- cumsum(h)
  ids <- stri_extract_first_regex(stri_sub(l[h], 2), "^\\S+")
  seqs <- vapply(split(l[!h], grp[!h]), paste, "", collapse = "")
  data.table(id = ids[as.integer(names(seqs))], seq = toupper(unname(seqs)))
}

write_fasta <- function(ids, seqs, path) {
  writeLines(paste0(">", ids, "\n", seqs), path)
}

## ---- Family assignment ------------------------------------------------------

family_of_class <- function(rclass) {
  fifelse(stri_detect_regex(rclass, "^SINE/Alu", case_insensitive = TRUE), "Alu",
  fifelse(stri_detect_regex(rclass, "^LINE/L1", case_insensitive = TRUE), "L1",
  fifelse(stri_detect_regex(rclass, "(^|/)SVA", case_insensitive = TRUE), "SVA", NA_character_)))
}

TAIL_REPEAT_RE <- "^\\((A|AA[A-Z]?|A{1,3}[CGT]?|T|TT[A-Z]?|T{1,3}[ACG]?)\\)n$|^[AT]-rich$"
SVA_HEXAMER_RE <- "^\\((CCCTCT|CCTCTC|CTCTCC|TCTCCC|CTCCCT|TCCCTC|AGAGGG|GAGGGA|AGGGAG|GGGAGA|GGAGAG|GAGAGG)\\)n$"

## ---- RepeatMasker line parsing ----------------------------------------------

## lines: character vector of RepeatMasker .out lines; ids: matching insertion id
parse_rm_lines <- function(lines, ids) {
  lines <- stri_trim_both(stri_replace_all_fixed(lines, "*", " "))
  tok <- stri_split_regex(lines, "\\s+")
  ok <- lengths(tok) >= 14
  tok <- tok[ok]; ids <- ids[ok]
  if (!length(tok)) {
    return(data.table(id = character(), qbeg = integer(), qend = integer(), strand = character(),
                      name = character(), rclass = character(), cbeg = integer(), cend = integer(),
                      cleft = integer(), div = numeric()))
  }
  m <- do.call(rbind, lapply(tok, `[`, 1:14))
  unparen <- function(x) suppressWarnings(as.integer(stri_replace_all_regex(x, "[()]", "")))
  plus <- m[, 9] == "+"
  dt <- data.table(
    id = ids, div = suppressWarnings(as.numeric(m[, 2])),
    qbeg = as.integer(m[, 6]), qend = as.integer(m[, 7]),
    strand = fifelse(plus, "+", "C"), name = m[, 10], rclass = m[, 11],
    c12 = unparen(m[, 12]), c13 = unparen(m[, 13]), c14 = unparen(m[, 14])
  )
  dt[, `:=`(cbeg = fifelse(plus, c12, c14), cend = fifelse(plus, c13, c13),
            cleft = fifelse(plus, c14, c12))]
  dt[, c("c12", "c13", "c14") := NULL]
  dt[, `:=`(lo = pmin(cbeg, cend), hi = pmax(cbeg, cend))][, `:=`(cbeg = lo, cend = hi)][, c("lo", "hi") := NULL]
  dt[!is.na(qbeg) & !is.na(qend) & qend >= qbeg & !is.na(cbeg)]
}

## Element_Hits cells -> long hits table
parse_hits_field <- function(cells, ids) {
  cells[is.na(cells) | cells %in% c("", "NONE", "nan", "NA")] <- ""
  lst <- stri_extract_all_regex(cells, "'[^']*'|\"[^\"]*\"")
  n <- lengths(lst)
  n[vapply(lst, function(x) length(x) == 1 && is.na(x[1]), TRUE)] <- 0L
  lines <- stri_sub(unlist(lst[n > 0]), 2, -2)
  parse_rm_lines(lines, rep(ids, n))
}

first_col <- function(cols, ...) {
  for (n in c(...)) {
    hit <- cols[tolower(cols) == tolower(n)]
    if (length(hit)) return(hit[1])
  }
  NULL
}

## L1ME-AID's ID is `chrom:POS;REF`, which two different ALT alleles at one
## site share, so every row gets its own key (`id`) and hits are keyed on it;
## the original ID is kept as `limeaid_id` for the allele-table join.
## L1ME-AID v1.4 column order, used when the table has no header (shards
## concatenated without one, as AnnotateL1MEAID's ConcatIntactMEIShards does).
## Older versions name two of these TE_Hits / FILTER_RESULTS instead.
LIMEAID_COLS <- c("ID", "Sequence", "Element_Hits", "Sequence_Length", "Element_Annotation",
                  "Element_Divergence", "Element_Proportion", "Element_Percentage", "Element_Designation",
                  "Orientation", "FLAGS", "Tail_Begins", "Tail_Type", "Tail_Length", "Tail_Seed_Hits",
                  "Unique_Element_Count", "Twin_Priming_Flag")

read_limeaid <- function(path) {
  d <- fread_any(path, sep = "\t", quote = "", header = FALSE, fill = TRUE, colClasses = "character",
                 showProgress = FALSE)
  has_header <- nrow(d) > 0 && any(tolower(unlist(d[1])) == "sequence")
  if (has_header) {
    setnames(d, make.unique(as.character(unlist(d[1]))))
    d <- d[-1]
  } else {
    n <- min(ncol(d), length(LIMEAID_COLS))
    setnames(d, seq_len(n), LIMEAID_COLS[seq_len(n)])
    message(path, ": no header; reading it as the L1ME-AID v1.4 column order (", ncol(d), " columns)")
  }
  layout <- paste0("expected L1ME-AID output, either with a header or headerless in the v1.4 order: ",
                   paste(LIMEAID_COLS, collapse = ", "))
  cols <- names(d)
  c_id <- first_col(cols, "ID"); if (is.null(c_id)) c_id <- cols[1]
  c_seq <- first_col(cols, "Sequence")
  c_hits <- first_col(cols, "Element_Hits", "TE_Hits")
  if (is.null(c_seq) || is.null(c_hits)) stop(path, ": no Sequence / Element_Hits columns; ", layout)
  ## Positions are only an assumption without a header, so check the content.
  s <- head(d[[c_seq]], 2000)
  h <- head(d[[c_hits]], 2000)
  seq_ok <- mean(stri_detect_regex(s, "^[ACGTNacgtn]+$"))
  hits_ok <- mean(stri_detect_regex(h, "^\\s*\\[") | h %in% c("", "NONE", "nan", "NA", "[]"))
  if (seq_ok < 0.95 || hits_ok < 0.95) {
    stop(sprintf("%s: column '%s' is %.0f%% DNA and column '%s' is %.0f%% RepeatMasker hit lists (need >= 95%% each); %s",
                 path, c_seq, 100 * seq_ok, c_hits, 100 * hits_ok, layout))
  }
  c_ann <- first_col(cols, "Element_Annotation")
  c_flag <- first_col(cols, "FLAGS", "FILTER_RESULTS")
  flag <- if (!is.null(c_flag)) d[[c_flag]] else NA_character_
  intact <- stri_detect_regex(flag, "^(INTACT|soloLTR|ERV)")
  key <- paste0("row", seq_len(nrow(d)))
  ins <- data.table(
    id = key, limeaid_id = d[[c_id]], seq = toupper(d[[c_seq]]),
    subfamily = if (!is.null(c_ann)) d[[c_ann]] else NA_character_,
    limeaid_flags = fifelse(intact %in% TRUE, NA_character_, flag),
    intact_flag = fifelse(intact %in% TRUE, flag, NA_character_)
  )
  hits <- parse_hits_field(d[[c_hits]], key)
  list(ins = ins, hits = hits)
}

## One table holding every family; the family is assigned after the
## allele-table join (assign_family).
load_insertions <- function(path) {
  d <- read_limeaid(path)
  d$hits[, fam_class := family_of_class(rclass)]
  d$ins[, L := nchar(seq)]
  d
}

family_of_type <- function(x) {
  fifelse(stri_detect_regex(x, "alu", case_insensitive = TRUE), "Alu",
  fifelse(stri_detect_regex(x, "line|l1", case_insensitive = TRUE), "L1",
  fifelse(stri_detect_regex(x, "sva", case_insensitive = TRUE), "SVA", NA_character_)))
}

## Family from the allele table's MEI_type, else the RepeatMasker family with
## the most bp in the sequence; sequences with neither are dropped. A blank
## L1ME-AID subfamily is filled from the longest same-family hit.
assign_family <- function(ins, hits) {
  ins[, family := family_of_type(MEI_type)]
  n_type <- ins[!is.na(family), .N]
  dom <- hits[!is.na(fam_class), .(bp = sum(qend - qbeg + 1L)), by = .(id, fam_class)][order(id, -bp)][
    , .SD[1], by = id]
  ins[dom, on = "id", family := fifelse(is.na(family), i.fam_class, family)]
  n_none <- ins[is.na(family), .N]
  message(sprintf("family: %d from MEI_type, %d from RepeatMasker hits, %d with neither (dropped).",
                  n_type, nrow(ins) - n_type - n_none, n_none))
  ins <- ins[!is.na(family)]
  top <- hits[ins[, .(id, family)], on = "id", nomatch = 0][fam_class == family][
    order(-(qend - qbeg)), .SD[1], by = id][, .(id, top = name)]
  ins[top, on = "id", subfamily := fifelse(is.na(subfamily) | subfamily %in% c("", "NONE", "nan"), i.top, subfamily)]
  ins
}

## Clustered MEI allele table (the `*.clustered.bed` lineage: one row per
## distinct ALT allele at a site, `cluster_id`/`var_id` grouping alleles into
## one insertion). Joined to L1ME-AID rows on chrom + POS + md5(ALT), because
## L1ME-AID IDs are `chrom:POS;REF` and its Sequence is the VCF ALT; md5 alone
## collides across sites (identical full-length alleles recur genome-wide).
MEI_TABLE_COLS <- c("chrom", "start", "end", "REF", "ALT", "MEI_type", "MEI_subfamily", "length",
                    "orientation", "intact_status", "cluster_id", "var_id", "ALT_md5")

read_mei_table <- function(path) {
  d <- fread_any(path, sep = "\t", colClasses = "character", showProgress = FALSE)
  if (!all(c("chrom", "start") %in% names(d))) {
    if (ncol(d) != length(MEI_TABLE_COLS)) {
      stop(path, ": no header and ", ncol(d), " columns; expected a header or the ",
           length(MEI_TABLE_COLS), " columns ", paste(MEI_TABLE_COLS, collapse = ", "))
    }
    d <- fread_any(path, sep = "\t", colClasses = "character", header = FALSE, showProgress = FALSE)
    setnames(d, MEI_TABLE_COLS)
  }
  d
}

attach_mei_table <- function(ins, path, variant_col = "var_id") {
  m <- read_mei_table(path)
  if (!variant_col %in% names(m)) stop("mei_table has no '", variant_col, "' column (mei_variant_col)")
  m[, md5 := if ("ALT" %in% names(m)) as.character(openssl::md5(toupper(ALT))) else tolower(ALT_md5)]
  ## The table's own ALT_md5 (md5 of the ALT as written) is what genotype tables are keyed on
  m[, ALT_md5 := if ("ALT_md5" %in% names(m)) tolower(ALT_md5) else md5]
  m[, pos := as.integer(start)]
  if (!"MEI_type" %in% names(m)) m[, MEI_type := NA_character_]
  m[, variant_id := if (variant_col == "cluster_id") paste(MEI_type, cluster_id, sep = "_cl") else get(variant_col)]
  key <- unique(m[, .(chrom, pos, md5, variant_id, ALT_md5, MEI_type)])
  n_conflict <- key[, uniqueN(variant_id), by = .(chrom, pos, md5)][V1 > 1, .N]
  key <- key[, .SD[1], by = .(chrom, pos, md5)]

  loc <- stri_match_first_regex(ins$limeaid_id, "^(.+?):(\\d+)(?:;|$)")
  ins[, `:=`(chrom = loc[, 2], id_pos = as.integer(loc[, 3]), md5 = as.character(openssl::md5(seq)))]
  if (all(is.na(ins$id_pos))) stop("L1ME-AID IDs are not chrom:POS[;REF] (e.g. ", ins$limeaid_id[1], "); cannot join mei_table")
  ## Table start is VCF POS in cluster_mei_gt.R; also try 0-based BED starts.
  hits_by_off <- vapply(c(0L, -1L), function(o)
    nrow(key[ins[, .(chrom, pos = id_pos + o, md5)], on = .(chrom, pos, md5), nomatch = 0]), 1L)
  off <- c(0L, -1L)[which.max(hits_by_off)]
  ins[, pos := id_pos + off]
  ins <- key[ins, on = .(chrom, pos, md5)]

  funnel <- ins[, .(l1meaid_rows = .N, matched = sum(!is.na(variant_id)),
                    variants = uniqueN(variant_id[!is.na(variant_id)]))]
  cat(sprintf("mei_table: %d allele rows, %d variants; joined with start = POS%s.\n",
              nrow(m), uniqueN(m$variant_id), if (off) " - 1" else ""))
  if (n_conflict) cat(sprintf("mei_table: %d alleles map to more than one variant; kept the first.\n", n_conflict))
  print(knitr::kable(funnel, caption = "L1ME-AID rows matched to mei_table (unmatched rows are dropped)"))
  if (!any(!is.na(ins$variant_id))) stop("No L1ME-AID row matched mei_table on chrom + POS + md5(ALT)")
  ins <- ins[!is.na(variant_id)]
  ins[, `:=`(id_pos = NULL, md5 = NULL)]
  ins
}

## Per-individual genotype table (CHROM, POS, ALT_md5, indiv, GT, hap1, hap2,
## ...), one row per (allele, individual), keyed like the clustered allele
## table (POS = its `start`). hap1/hap2 say which haplotypes carry the allele;
## GT is used only where they are missing.
read_genotypes <- function(path) {
  g <- fread_any(path, sep = "\t", colClasses = "character", showProgress = FALSE)
  need <- c("CHROM", "POS", "ALT_md5", "indiv")
  if (!all(need %in% names(g))) stop(path, ": genotype table needs columns ", paste(need, collapse = ", "))
  as_hap <- function(x) {
    x <- toupper(stri_trim_both(x))
    fifelse(x %in% c("TRUE", "T", "1"), TRUE, fifelse(x %in% c("FALSE", "F", "0"), FALSE, NA))
  }
  h1 <- if ("hap1" %in% names(g)) as_hap(g$hap1) else rep(NA, nrow(g))
  h2 <- if ("hap2" %in% names(g)) as_hap(g$hap2) else rep(NA, nrow(g))
  gt <- if ("GT" %in% names(g)) stri_split_regex(g$GT, "[/|]", simplify = TRUE) else matrix(NA_character_, nrow(g), 2)
  h1[is.na(h1)] <- (gt[is.na(h1), 1] %in% c("1"))
  h2[is.na(h2)] <- (if (ncol(gt) > 1) gt[is.na(h2), 2] else NA_character_) %in% c("1")
  d <- data.table(chrom = g$CHROM, pos = as.integer(g$POS), ALT_md5 = tolower(g$ALT_md5), indiv = g$indiv,
                  hap1 = h1, hap2 = h2)
  ## a row repeated for the same allele and individual counts once
  d[, .(hap1 = any(hap1), hap2 = any(hap2)), by = .(chrom, pos, ALT_md5, indiv)]
}

wmedian <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x)
  as.numeric(x[o][which(cumsum(w[o]) >= sum(w) / 2)[1]])
}

## ---- Orientation ------------------------------------------------------------

## Re-express each insertion's family hits in element-sense coordinates (the
## strand carrying most family bp becomes '+'); returns list(ins with
## sense/strand columns, fam_hits, other_hits), all in sense coordinates.
to_sense <- function(ins, hits) {
  h <- hits[ins[, .(id, family, L)], on = "id", nomatch = 0]
  h[, is_core := !is.na(fam_class) & fam_class == family]
  h[, is_fam := is_core | (family == "SVA" & stri_detect_regex(name, SVA_HEXAMER_RE, case_insensitive = TRUE))]
  cov <- h[is_core == TRUE, .(bp = sum(qend - qbeg + 1)), by = .(id, strand)][order(id, -bp)][, .SD[1], by = id]
  ins[, strand := "+"]
  ins[cov, on = "id", strand := fifelse(i.strand == "+", "+", "-")]
  minus <- ins[strand == "-", id]
  ins[, sense := seq]
  ins[strand == "-", sense := revcomp(seq)]
  h[id %in% minus, `:=`(qbeg = L - qend + 1L, qend = L - qbeg + 1L, strand = fifelse(strand == "+", "C", "+"))]
  h[, qlen := qend - qbeg + 1L]
  list(ins = ins, fam = h[is_fam == TRUE], other = h[is_fam == FALSE])
}

ARGS <- list(
  limeaid_table = c("", "Required. L1ME-AID/INTACT_MEI TSV covering every family"),
  mei_table = c("", "Required. Clustered MEI allele table (chrom, start, ..., ALT, MEI_type, ..., cluster_id, var_id, ALT_md5); groups sequences into variants and sets their family"),
  mei_variant_col = c("var_id", "mei-table column that defines a variant"),
  genotype_table = c("", "Required. Per-individual genotypes (CHROM, POS, ALT_md5, indiv, GT, hap1, hap2, ...), keyed like --mei-table"),
  n_individuals = c("0", "Cohort size for AF (AN = 2 x this); 0 = distinct indiv in --genotype-table"),
  max_5p_missing = c("Alu=20,L1=50,SVA=50", "Consensus bp missing at the 5' end before calling 5' truncation, per family"),
  max_3p_missing = c("Alu=20,L1=50,SVA=50", "Consensus bp missing at the 3' end before calling 3' truncation, per family"),
  min_transduction = c("30", "Minimum unique bp for a transduction"),
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
  transduction_paf = c("", "Precomputed PAF of transductions.fa vs the reference (instead of --reference)"),
  ref_rmsk = c("", "Reference repeats: UCSC rmsk.txt[.gz], RepeatMasker .out, or BED"),
  known_sources = c("", "Comma-separated BEDs of extra (e.g. non-reference) source elements"),
  genes = c("", "Gene/exon BED to annotate source loci"),
  minimap2_preset = c("sr", "minimap2 -x preset"),
  threads = c("4", "minimap2 threads"),
  min_mapq = c("20", "MAPQ for a unique source mapping"),
  min_identity = c("0.9", "Minimum identity of a source hit"),
  min_qcov = c("0.7", "Minimum query coverage of a source hit"),
  self_window = c("1000", "Hits this close to the insertion site are self hits"),
  source_window = c("2000", "Maximum distance from a source element"),
  overlap_tolerance = c("50", "Allowed overlap with a source element"),
  cluster_window = c("5000", "Window for grouping element-less source loci"),
  full_length_l1_bp = c("5900", "L1 source length counted as full-length"),
  run_ava = c("true", "All-vs-all cluster unmapped segments into shared non-reference sources"),
  outdir = c(".", "Output directory"),
  prefix = c("mei", "Output file prefix")
)
ap <- arg_parser("Catalog MEI structural forms and trace transduction sources")
for (n in names(ARGS)) ap <- add_argument(ap, paste0("--", gsub("_", "-", n)), help = ARGS[[n]][2], default = ARGS[[n]][1])
argv <- parse_args(ap)
params <- setNames(lapply(names(ARGS), function(n) {
  v <- argv[[n]]
  if (is.null(v) || (length(v) == 1 && is.na(v))) "" else as.character(v)
}), names(ARGS))

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
tab_path <- function(name, ext = "tsv") file.path(params$outdir, paste0(params$prefix, ".structure_", name, ".", ext))
save_tab <- function(d, name, ext = "tsv") {
  f <- tab_path(name, ext)
  fwrite(d, f, sep = "\t", na = "NA")
  message("wrote ", f, " (", nrow(d), " rows)")
}
LOG <- function(...) message(sprintf(...))

P <- list(
  max5 = lapply(p_map(params$max_5p_missing), as.integer),
  max3 = lapply(p_map(params$max_3p_missing), as.integer),
  min_td = p_int(params$min_transduction), max_td_masked = p_num(params$max_transduction_masked),
  min_orphan = p_int(params$min_orphan), min_polya = p_int(params$min_polya),
  polya_frac = p_num(params$polya_min_frac), polya_gap = p_int(params$polya_gap),
  slack = p_int(params$terminal_slack), min_inv = p_int(params$min_inversion),
  inv_slack = p_int(params$inv_slack), min_gap = p_int(params$min_internal_gap)
)
for (req in c("limeaid_table", "mei_table", "genotype_table")) {
  if (is.null(p_str(params[[req]]))) stop("--", gsub("_", "-", req), " is required")
}
run_sources <- !is.null(p_str(params$reference)) || !is.null(p_str(params$transduction_paf))

## ---- Load -------------------------------------------------------------------

dat <- load_insertions(params$limeaid_table)
hits <- dat$hits
ins <- attach_mei_table(dat$ins, params$mei_table, params$mei_variant_col)
ins <- assign_family(ins, hits)
for (f in unique(ins$family)) {
  if (is.null(P$max5[[f]])) P$max5[[f]] <- 50L
  if (is.null(P$max3[[f]])) P$max3[[f]] <- 50L
}
## Identical (family, sequence, hits) inputs give identical calls: classify
## each distinct one once. Hit signature excludes the query name.
sig <- hits[, .(sig = paste(qbeg, qend, strand, name, cbeg, cend, collapse = ";")), by = id]
ins[sig, on = "id", sig := i.sig]
ins[is.na(sig), sig := ""]
ins[, uid := .GRP, by = .(family, seq, sig)]
rep_ids <- ins[, .(id = id[1]), by = uid]
LOG("Loaded %d sequences (%d distinct) across %s.", nrow(ins), nrow(rep_ids),
    paste(sprintf("%s: %d", names(table(ins$family)), table(ins$family)), collapse = ", "))

## ---- Classification ---------------------------------------------------------

SNRNA_RE <- "snRNA|scRNA|srpRNA|7SL|U6|U3|7SK"
SPLIT_3P_RE <- "_(3end|orf2)$"
PAS_RE <- "AATAAA|ATTAAA|AGTAAA|TATAAA|CATAAA|GATAAA|AATATA|AATACA|AATAGA|ACTAAA|AAGAAA|AATGAA"

## Fraction of each interval [s, e) (0-based) covered by hits (qbeg..qend, 1-based).
masked_frac <- function(iv, h) {
  if (!nrow(iv)) return(numeric())
  iv <- copy(iv)[, row := .I]
  if (!nrow(h)) return(rep(0, nrow(iv)))
  ov <- h[, .(id, qbeg, qend)][iv, on = "id", nomatch = 0, allow.cartesian = TRUE][
    , .(ov = sum(pmax(0L, pmin(e, qend) - pmax(s, qbeg - 1L)))), by = row]
  out <- numeric(nrow(iv))
  out[ov$row] <- ov$ov
  pmin(1, out / pmax(1, iv$e - iv$s))
}

## Maximal A-rich segments (0-based half-open, absolute) in sense[ws+1 .. end]:
## runs of A merged across <= 2 non-A bases; kept if >= min_len and >= frac A.
## Leading sparse runs are trimmed off so they cannot sink a real tail.
find_tails <- function(ids, seqs, ws, min_len, frac) {
  win <- stri_sub(seqs, ws + 1L)
  loc <- stri_locate_all_regex(win, "A+(?:[^A]{1,2}A+)*", omit_no_match = TRUE)
  n <- vapply(loc, nrow, 1L)
  if (!sum(n)) return(data.table(id = character(), s = integer(), e = integer()))
  m <- do.call(rbind, loc)
  b <- data.table(id = rep(ids, n), ws = rep(ws, n), st = m[, 1], en = m[, 2])
  b[, str := stri_sub(rep(win, n), st, en)]
  b <- b[nchar(str) >= min_len]
  ok <- function(x) nchar(x) >= min_len & stri_count_fixed(x, "A") / nchar(x) >= frac
  b[, pass := ok(str)]
  for (it in 1:100) {
    todo <- which(!b$pass & stri_detect_regex(b$str, "^A+[^A]{1,2}A"))
    if (!length(todo)) break
    trimmed <- stri_replace_first_regex(b$str[todo], "^A+[^A]{1,2}", "")
    b[todo, `:=`(st = st + (nchar(str) - nchar(trimmed)), str = trimmed)]
    b[todo, pass := ok(str)]
  }
  b[pass == TRUE, .(id, s = ws + st - 1L, e = ws + en)]
}

classify_structures <- function(ins, hits, P) {
  sx <- to_sense(ins, hits)
  I <- sx$ins
  core <- sx$fam[is_core == TRUE]
  main <- core[strand == "+"]
  inv <- core[strand == "C"]
  hex <- sx$fam[is_core == FALSE]
  is_tail <- stri_detect_regex(sx$other$name, TAIL_REPEAT_RE, case_insensitive = TRUE)
  tailh <- sx$other[is_tail & strand == "+"]
  nontail <- sx$other[!is_tail]

  R <- I[, .(id, family, subfamily, L, strand, limeaid_flags, intact_flag, sense)]
  R <- merge(R, core[, .(family_bp = sum(qlen)), by = id], by = "id", all.x = TRUE)
  R[is.na(family_bp), family_bp := 0L][, family_frac := family_bp / L]

  setorder(main, id, qbeg)
  agg <- main[, .(sense_start = min(qbeg) - 1L, te_end = max(qend), retained = sum(cend - cbeg + 1L),
                  divergence = sum(div * qlen) / sum(qlen)), by = id]
  five <- main[main[, .I[which.min(qbeg)], by = id]$V1, .(id, cons_start = cbeg, five_hit = name, five_cend = cend)]
  three <- main[main[, .I[which.max(qend)], by = id]$V1,
                .(id, cons_end = cend, cons_left_3p = cleft, cons_len = cend + cleft, three_hit = name)]
  tst <- rbind(main[, .(id, qbeg)], hex[, .(id, qbeg)], inv[, .(id, qbeg)])[, .(te_start = min(qbeg) - 1L), by = id]
  R <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(R, agg, five, three, tst))
  R[, has_main := !is.na(te_end)]

  ## b/c/i: truncation
  R[, max5 := unlist(P$max5[family])][, max3 := unlist(P$max3[family])]
  R[, five_prime_missing := cons_start - 1L]
  R[, five_prime_truncated := has_main & (stri_detect_regex(five_hit, SPLIT_3P_RE, case_insensitive = TRUE) |
                                           five_prime_missing > max5)]
  R[, retained_cons_frac := pmin(1, retained / cons_len)]
  R[, three_prime_missing := cons_left_3p]
  R[, three_prime_truncated := has_main & cons_left_3p > max3 & !stri_detect_regex(three_hit, "_5end$")]

  ## e: 5' inversion (opposite-strand family sequence 5' of the sense portion)
  inv5 <- inv[R[has_main == TRUE, .(id, sense_start)], on = "id", nomatch = 0][qend <= sense_start + P$inv_slack]
  if (nrow(inv5)) {
    ia <- inv5[, .(inversion_len = sum(qlen), inversion_cons_start = min(cbeg), inversion_cons_end = max(cend),
                   inner_cend = cend[which.max(qend)]), by = id]
    R <- merge(R, ia, by = "id", all.x = TRUE)
  } else {
    R[, `:=`(inversion_len = NA_integer_, inversion_cons_start = NA_integer_,
             inversion_cons_end = NA_integer_, inner_cend = NA_integer_)]
  }
  R[is.na(inversion_len), inversion_len := 0L]
  R[, five_prime_inversion := inversion_len >= P$min_inv]
  ## junction > 0: deletion at the inversion point, < 0: duplication (Ostertag & Kazazian 2001)
  R[five_prime_inversion == TRUE, inversion_junction := cons_start - inner_cend - 1L]
  if (!"inversion_junction" %in% names(R)) R[, inversion_junction := NA_integer_]
  R[, inner_cend := NULL]

  ## internal deletions / rearrangements between consecutive sense hits
  main[, `:=`(nb = shift(cbeg, -1), nq = shift(qbeg, -1), nn = shift(name, -1)), by = id]
  pairs <- main[!is.na(nb) & (name == nn | stri_detect_regex(nn, SPLIT_3P_RE) | stri_detect_regex(name, "_5end$"))]
  pairs[, `:=`(cgap = nb - cend - 1L, qgap = nq - qend - 1L)]
  pd <- pairs[, .(internal_deletion_bp = sum(pmax(0L, (cgap - pmax(qgap, 0L))) * ((cgap - pmax(qgap, 0L)) >= P$min_gap)),
                  internal_rearrangement = any(cgap < -P$min_gap)), by = id]
  R <- merge(R, pd, by = "id", all.x = TRUE)
  R[is.na(internal_deletion_bp), internal_deletion_bp := 0L]
  R[, internal_deletion := internal_deletion_bp > 0]
  R[is.na(internal_rearrangement), internal_rearrangement := FALSE]

  ## snRNA chimeras, other repeats
  sn <- nontail[stri_detect_regex(rclass, SNRNA_RE) | stri_detect_regex(name, SNRNA_RE)]
  R <- merge(R, sn[, .(snRNA_name = paste(sort(unique(name)), collapse = ",")), by = id], by = "id", all.x = TRUE)
  R[, snRNA_chimera := !is.na(snRNA_name) & has_main]
  ot <- nontail[!(id %in% sn$id & name %in% sn$name) & !rclass %in% c("Simple_repeat", "Low_complexity")]
  R <- merge(R, ot[, .(other_repeat_bp = sum(qlen)), by = id], by = "id", all.x = TRUE)
  R[is.na(other_repeat_bp), other_repeat_bp := 0L]

  ## poly(A) tails: search from just before the element's 3' end (or the
  ## whole sequence when there is no sense element)
  R[, ws := fifelse(has_main, pmax(0L, te_end - 40L), 0L)]
  tails <- rbind(find_tails(R$id, R$sense, R$ws, P$min_polya, P$polya_frac),
                 tailh[, .(id, s = qbeg - 1L, e = qend)])
  tails <- unique(tails)[R[, .(id, L, te_end, ws)], on = "id", nomatch = 0][e > ws]
  term <- tails[e >= L - P$slack][order(id, -(e - s))][, .SD[1], by = id][, .(id, term_s = s, term_e = e)]
  first <- tails[!is.na(te_end) & s <= te_end + P$polya_gap & e > te_end - 10L][order(id, s)][, .SD[1], by = id][
    , .(id, first_s = s, first_e = e)]
  after <- tails[!is.na(te_end) & s <= te_end + P$polya_gap & e >= te_end][order(id, s)][, .SD[1], by = id][
    , .(id, after_s = s)]
  R <- Reduce(function(a, b) merge(a, b, by = "id", all.x = TRUE), list(R, term, first, after))
  R[, terminal_polyA := !is.na(term_s)]
  R[, polyA_len := fifelse(terminal_polyA, term_e - term_s, 0L)]
  R[, internal_polyA := has_main & !is.na(first_s) & terminal_polyA & first_e < term_s]

  ## d: 3' transduction
  R[, `:=`(tdx3_start = NA_integer_, tdx3_end = NA_integer_, tdx3_masked_frac = NA_real_, tdx3_type = NA_character_)]
  R[has_main & terminal_polyA, `:=`(tdx3_start = fifelse(internal_polyA, first_e, te_end), tdx3_end = term_s)]
  R[, tdx3_len := tdx3_end - tdx3_start]
  cand3 <- R[!is.na(tdx3_len) & tdx3_len >= P$min_td, .(id, s = tdx3_start, e = tdx3_end)]
  if (nrow(cand3)) R[cand3, on = "id", tdx3_masked_frac := masked_frac(cand3, nontail)]
  R[, three_prime_transduction := !is.na(tdx3_masked_frac) & tdx3_masked_frac <= P$max_td_masked]
  R[three_prime_transduction == FALSE, `:=`(tdx3_start = NA_integer_, tdx3_end = NA_integer_, tdx3_len = NA_integer_)]
  R[three_prime_transduction == TRUE, tdx3_type := fifelse(internal_polyA, "TE-A-unique-A", "TE-unique-A")]
  R[, unresolved_3p_extension := has_main & !terminal_polyA & (L - te_end - P$slack) >= P$min_td]

  ## orphan candidates: no sense element, unique sequence + terminal poly(A)
  R[, orphan_candidate := FALSE]
  orph <- R[has_main == FALSE & terminal_polyA, .(id, s = 0L, e = term_s)]
  if (nrow(orph)) {
    orph[, uniq := (e - s) * (1 - masked_frac(orph, sx$other))]
    R[orph, on = "id", orphan_candidate := i.uniq >= P$min_orphan]
    R[orphan_candidate == TRUE, `:=`(tdx3_start = 0L, tdx3_end = term_s, tdx3_len = term_s)]
  }

  ## f/h: 5' foreign sequence -> 5' transduction (intact element 5' end) or
  ## exon-trapping candidate (internally starting element, spliced junction)
  R[, `:=`(five_prime_transduction = FALSE, exon_trapping_candidate = FALSE,
           tdx5_len = NA_integer_, tdx5_end = NA_integer_, tdx5_masked_frac = NA_real_,
           tdx5_t_rich_head = NA, tdx5_junction_dinuc = NA_character_)]
  cand5 <- R[has_main & (te_start - P$slack) >= P$min_td, .(id, s = P$slack, e = te_start)]
  if (nrow(cand5)) {
    R[cand5, on = "id", tdx5_masked_frac := masked_frac(cand5, nontail)]
    R[id %in% cand5$id, tdx5_t_rich_head := stri_count_fixed(stri_sub(sense, 1, 20), "T") / 20 >= 0.8]
    ok5 <- R$id %in% cand5$id & R$tdx5_masked_frac <= P$max_td_masked & !R$tdx5_t_rich_head
    ok5[is.na(ok5)] <- FALSE
    R[ok5 & five_prime_truncated & !five_prime_inversion, exon_trapping_candidate := TRUE]
    R[ok5 & !exon_trapping_candidate, five_prime_transduction := TRUE]
    R[exon_trapping_candidate == TRUE, tdx5_junction_dinuc := stri_sub(sense, pmax(1L, te_start - 1L), te_start)]
    R[five_prime_transduction | exon_trapping_candidate, `:=`(tdx5_len = te_start, tdx5_end = te_start)]
  }

  ## g/i: mechanism of 3' truncation
  R[, `:=`(three_prime_end_A_frac = NA_real_, pas_upstream_of_tail = NA)]
  R[three_prime_truncated == TRUE, three_prime_end_A_frac :=
      stri_count_fixed(stri_sub(sense, pmax(1L, te_end - 19L), te_end), "A") /
      nchar(stri_sub(sense, pmax(1L, te_end - 19L), te_end))]
  R[three_prime_truncated == TRUE, pas_upstream_of_tail := {
    anchor <- fifelse(is.na(after_s), te_end, after_s)
    stri_detect_regex(stri_sub(sense, pmax(1L, anchor - 39L), pmax(0L, anchor - 5L)), PAS_RE)
  }]
  R[, internal_priming := three_prime_truncated & !is.na(three_prime_end_A_frac) & three_prime_end_A_frac >= 0.7]
  R[, premature_polyadenylation := three_prime_truncated & !internal_priming & !is.na(after_s) & pas_upstream_of_tail %in% TRUE]

  ## combined label and transduced sequences
  R[, full_length := has_main & !five_prime_truncated & !five_prime_inversion & !three_prime_truncated]
  flag_order <- c("five_prime_truncated", "five_prime_inversion", "three_prime_truncated", "internal_priming",
                  "premature_polyadenylation", "internal_deletion", "internal_rearrangement",
                  "three_prime_transduction", "five_prime_transduction", "exon_trapping_candidate", "snRNA_chimera")
  lab <- do.call(cbind, lapply(flag_order, function(f) fifelse(R[[f]] %in% TRUE, f, NA_character_)))
  lab <- cbind(lab, fifelse(R$has_main & !R$terminal_polyA, "no_terminal_polyA", NA_character_))
  R[, structure := apply(lab, 1, function(z) paste(z[!is.na(z)], collapse = "+"))]
  R[structure == "", structure := "full_length_canonical"]
  R[has_main == FALSE, structure := "no_element"]
  R[, `:=`(tdx3_seq = NA_character_, tdx5_seq = NA_character_)]
  R[!is.na(tdx3_len), tdx3_seq := stri_sub(sense, tdx3_start + 1L, tdx3_end)]
  R[five_prime_transduction | exon_trapping_candidate, tdx5_seq := stri_sub(sense, 1L, te_start)]
  R[, c("max5", "max3", "ws", "term_s", "term_e", "first_s", "first_e", "after_s", "sense", "five_cend") := NULL]
  R
}

FLAGS <- c("full_length", "five_prime_truncated", "three_prime_transduction", "five_prime_inversion",
           "five_prime_transduction", "internal_priming", "exon_trapping_candidate", "three_prime_truncated",
           "premature_polyadenylation", "internal_deletion", "internal_rearrangement", "orphan_candidate",
           "snRNA_chimera", "terminal_polyA", "internal_polyA")

res <- classify_structures(ins[.(rep_ids$id), on = "id"], hits[id %in% rep_ids$id], P)
res <- merge(rep_ids, res, by = "id")[, id := NULL]
per_seq <- merge(ins[, .(id, limeaid_id, variant_id, chrom, pos, ALT_md5, uid)][, seq_md5 := uid],
                 res, by = "uid")
per_seq[, uid := NULL]
for (f in FLAGS) per_seq[is.na(get(f)), (f) := FALSE]

## Genotypes: haplotypes carrying each catalogued allele. A variant's AF is
## its alleles' carried haplotypes over AN = 2 x cohort size.
G <- read_genotypes(params$genotype_table)
n_ind <- p_int(params$n_individuals)
if (is.na(n_ind) || n_ind <= 0) n_ind <- uniqueN(G$indiv)
AN <- 2L * n_ind
G <- G[hap1 | hap2]
ac <- G[, .(allele_count = sum(hap1) + sum(hap2), n_carriers = .N), by = .(chrom, pos, ALT_md5)]
per_seq[ac, on = .(chrom, pos, ALT_md5), `:=`(allele_count = i.allele_count, n_carriers = i.n_carriers)]
per_seq[is.na(allele_count), `:=`(allele_count = 0L, n_carriers = 0L)]
per_seq[, af := allele_count / AN]
carried <- G[unique(per_seq[, .(chrom, pos, ALT_md5, variant_id)]), on = .(chrom, pos, ALT_md5), nomatch = 0]
var_gt <- carried[, .(n_carriers = uniqueN(indiv), allele_count = sum(hap1) + sum(hap2)), by = variant_id]
n_unmatched <- nrow(ac[!per_seq, on = .(chrom, pos, ALT_md5)])
LOG("genotype table: %d individuals (AN = %d); %d of %d catalogued alleles carried by someone; %d genotyped alleles not in the catalog.",
    n_ind, AN, per_seq[allele_count > 0, .N], nrow(per_seq), n_unmatched)
if (!per_seq[allele_count > 0, .N]) stop("No catalogued allele matched the genotype table on CHROM + POS + ALT_md5")
save_tab(per_seq[, !c("tdx3_seq", "tdx5_seq"), with = FALSE], "per_sequence", "tsv.gz")

## Variant-level rollup: modal structure, majority flags, medians and
## concordance, each sequence weighted by the haplotypes carrying it (a
## variant none of whose alleles is carried falls back to equal weights).
per_seq[, w := as.numeric(allele_count)]
per_seq[, w := if (sum(w) > 0) w else rep(1, .N), by = variant_id]
wmode <- function(x, w) { s <- tapply(w, x, sum); names(s)[which.max(s)] }
num_cols <- c("L", "five_prime_missing", "three_prime_missing", "retained_cons_frac", "polyA_len",
              "tdx3_len", "tdx5_len", "inversion_len", "divergence")
per_var <- per_seq[, c(
  list(family = family[1], subfamily = wmode(fifelse(is.na(subfamily), "NA", subfamily), w),
       chrom = chrom[1], pos = pos[1], n_sequences = .N, n_distinct_seqs = uniqueN(seq_md5),
       structure = wmode(structure, w)),
  lapply(.SD[, FLAGS, with = FALSE], function(x) sum(w * x) / sum(w) >= 0.5),
  lapply(.SD[, num_cols, with = FALSE], function(x) wmedian(x, w))
), by = variant_id]
conc <- per_seq[per_var[, .(variant_id, modal = structure)], on = "variant_id"][
  , .(structure_concordance = sum(w * (structure == modal)) / sum(w), allele_concordance = mean(structure == modal)),
  by = variant_id]
per_var <- merge(per_var, conc, by = "variant_id")
per_var <- merge(per_var, var_gt, by = "variant_id", all.x = TRUE)
per_var[is.na(n_carriers), `:=`(n_carriers = 0L, allele_count = 0L)]
per_var[, af := allele_count / AN]
per_seq[, w := NULL]
save_tab(per_var, "per_variant")
LOG("%d sequences -> %d variants.", nrow(per_seq), nrow(per_var))

## One transduced segment per variant (the one carried by most sequences),
## named variant|kind|family|chrom|pos for source tracing.
tdx <- rbind(
  per_seq[!is.na(tdx3_seq), .(kind = "3p", seq = tdx3_seq, family, chrom, pos), by = variant_id],
  per_seq[!is.na(tdx5_seq), .(kind = "5p", seq = tdx5_seq, family, chrom, pos), by = variant_id]
)[, .N, by = .(variant_id, kind, seq, family, chrom, pos)][order(-N)][, .SD[1], by = .(variant_id, kind)]
tdx[, name := paste(variant_id, kind, family, chrom, pos, sep = "|")]
tdx_fa <- tab_path("transductions", "fa")
write_fasta(tdx$name, tdx$seq, tdx_fa)
LOG("wrote %s (%d segments)", tdx_fa, nrow(tdx))

if (!run_sources) {
  LOG("No --reference / --transduction-paf: skipping transduction source tracing.")
  quit(save = "no", status = 0)
}
if (!nrow(tdx)) {
  LOG("No transduced segments: nothing to trace.")
  quit(save = "no", status = 0)
}

## ---- Transduction source tracing --------------------------------------------
## Segments are mapped to the reference; hits within self_window of the
## insertion's own site are flanking-sequence capture, not transduction. Each
## remaining hit is assigned to a same-family, same-strand element upstream
## (3') or downstream (5') of it within source_window. Unmapped/self-only
## segments that recur are clustered all-vs-all: shared non-reference sources.

read_paf <- function(path) {
  l <- readLines(path)
  if (!length(l)) return(data.table())
  t <- stri_split_fixed(l, "\t", n = 13, simplify = TRUE)
  d <- data.table(qname = t[, 1], qlen = as.integer(t[, 2]), qstart = as.integer(t[, 3]), qend = as.integer(t[, 4]),
                  strand = t[, 5], tname = t[, 6], tlen = as.integer(t[, 7]), tstart = as.integer(t[, 8]),
                  tend = as.integer(t[, 9]), matches = as.integer(t[, 10]), alen = as.integer(t[, 11]),
                  mapq = as.integer(t[, 12]))
  d[, primary := stri_detect_fixed(t[, 13], "tp:A:P")]
  d
}
read_elements <- function(path, set) {
  if (grepl("\\.out(\\.gz)?$", path)) {                       # RepeatMasker .out (whitespace-aligned)
    l <- if (grepl("\\.gz$", path)) readLines(gzfile(path)) else readLines(path)
    l <- l[stri_detect_regex(l, "^\\s*\\d")]
    t <- stri_split_regex(stri_trim_both(l), "\\s+", simplify = TRUE)
    e <- data.table(chrom = t[, 5], start = as.integer(t[, 6]) - 1L, end = as.integer(t[, 7]),
                    strand = fifelse(t[, 9] == "+", "+", "-"), name = t[, 10], rclass = t[, 11])
  } else {
    d <- fread_any(path, header = FALSE, fill = TRUE, colClasses = "character", showProgress = FALSE)
    if (ncol(d) >= 17 && all(d$V10 %in% c("+", "-"))) {         # UCSC rmsk.txt
      e <- d[, .(chrom = V6, start = as.integer(V7), end = as.integer(V8), strand = V10, name = V11,
                 rclass = paste0(V12, "/", V13))]
    } else {                                                     # BED
      e <- d[, .(chrom = V1, start = as.integer(V2), end = as.integer(V3),
                 name = if (ncol(d) > 3) V4 else ".", strand = if (ncol(d) > 5) V6 else ".",
                 rclass = if (ncol(d) > 6) V7 else "")]
    }
  }
  e[, set := set]
  e[!is.na(start)]
}

threads <- params$threads
tseqs <- read_fasta(tdx_fa)
if (!nrow(tseqs)) tseqs <- data.table(id = character(), seq = character())
paf <- p_str(params$transduction_paf)
if (is.null(paf) && nrow(tseqs)) {
  paf <- tab_path("transductions", "paf")
  mm2 <- which_tool("minimap2"); if (is.null(mm2)) stop("minimap2 not found on PATH")
  LOG("Mapping %d transduced segments with minimap2 -x %s", nrow(tseqs), params$minimap2_preset)
  run_tool(mm2, c("-c", "-x", params$minimap2_preset, "--secondary=yes", "-N", "5", "-t", threads,
                  p_str(params$reference), tdx_fa), stdout = paf)
}
H <- if (!is.null(paf) && file.exists(paf)) read_paf(paf) else data.table()
meta <- tseqs[, {s <- stri_split_fixed(id, "|", simplify = TRUE)
  .(transduction = id, variant_id = s[, 1], kind = s[, 2], family = s[, 3], ins_chrom = s[, 4],
    ins_pos = suppressWarnings(as.integer(s[, 5])), td_len = nchar(seq))}]
best <- data.table(qname = character(), tname = character(), tstart = integer(), tend = integer(),
                   strand = character(), mapq = integer(), identity = numeric(), qcov = numeric(), n_good_hits = integer())
meta[, maps_to_insertion_site := FALSE]
if (nrow(H)) {
  H[, `:=`(identity = matches / alen, qcov = (qend - qstart) / qlen)]
  H <- H[identity >= p_num(params$min_identity) & qcov >= p_num(params$min_qcov)]
  H <- merge(H, meta[, .(qname = transduction, ins_chrom, ins_pos)], by = "qname")
  sw <- p_int(params$self_window)
  H[, self_hit := tname == ins_chrom & !is.na(ins_pos) & (abs(tstart - ins_pos) < sw | abs(tend - ins_pos) < sw)]
  meta[, maps_to_insertion_site := transduction %in% H[self_hit == TRUE, qname]]
  H <- H[self_hit == FALSE]
  if (nrow(H)) {
    best <- H[order(qname, -primary, -mapq, -matches)][, .SD[1], by = qname]
    best[, n_good_hits := H[, .N, by = qname][best, on = "qname", N]]
  }
}
src <- merge(meta, best[, .(transduction = qname, src_chrom = tname, src_start = tstart, src_end = tend,
                            src_strand = strand, mapq, identity, qcov, n_good_hits)], by = "transduction", all.x = TRUE)
src[, mapping := fifelse(is.na(src_chrom), fifelse(maps_to_insertion_site, "self_only", "unmapped"),
                 fifelse(mapq >= p_int(params$min_mapq), "unique", "multi"))]

## Source element: 3' transduction -> element ends just UPSTREAM of the
## segment (in the segment's genomic orientation); 5' -> starts just DOWNSTREAM.
E <- rbindlist(c(
  if (!is.null(p_str(params$ref_rmsk))) list(read_elements(params$ref_rmsk, "reference")),
  lapply(Filter(nzchar, stri_split_fixed(params$known_sources, ",")[[1]]), function(b) read_elements(b, basename(b)))
), fill = TRUE)
src[, `:=`(source_element = NA_character_, source_class = NA_character_, source_set = NA_character_,
           source_chrom = NA_character_, source_start = NA_integer_, source_end = NA_integer_,
           source_strand = NA_character_, source_distance = NA_integer_, source_full_length = NA)]
win <- p_int(params$source_window); tol <- p_int(params$overlap_tolerance)
FAM_RE <- c(Alu = "^SINE/Alu|(^|/)Alu", L1 = "^LINE/L1|(^|/)L1", SVA = "(^|/)SVA")
if (nrow(E) && any(!is.na(src$src_chrom))) {
  E <- E[set != "reference" | stri_detect_regex(rclass, "Alu|L1|SVA")]
  m <- src[!is.na(src_chrom)]
  m[, upstream := (kind == "3p") == (src_strand == "+")]
  m[, need := fifelse(src_strand == "+", "+", "-")]
  m[, `:=`(lo = fifelse(upstream, src_start - win, src_end - tol), hi = fifelse(upstream, src_start + tol, src_end + win))]
  Eend <- E[, .(chrom, pos = end, el_start = start, el_end = end, el_strand = strand, el_name = name, el_class = rclass, set)]
  Est <- E[, .(chrom, pos = start, el_start = start, el_end = end, el_strand = strand, el_name = name, el_class = rclass, set)]
  j1 <- Eend[m[upstream == TRUE], on = .(chrom == src_chrom, pos >= lo, pos <= hi), nomatch = 0, allow.cartesian = TRUE,
             .(transduction, family, need, src_start, el_start, el_end, el_strand, el_name, el_class, set, dist = src_start - el_end)]
  j2 <- Est[m[upstream == FALSE], on = .(chrom == src_chrom, pos >= lo, pos <= hi), nomatch = 0, allow.cartesian = TRUE,
            .(transduction, family, need, src_start, el_start, el_end, el_strand, el_name, el_class, set, dist = el_start - src_end)]
  j <- rbind(j1, j2)
  j <- j[el_strand %in% c(need, ".") & dist >= -tol & dist <= win]
  j <- j[set != "reference" | stri_detect_regex(paste0(el_class, "/", el_name), FAM_RE[family])]
  j <- j[order(transduction, abs(dist), -(el_end - el_start))][, .SD[1], by = transduction]
  src[j, on = "transduction", `:=`(source_element = i.el_name, source_class = i.el_class, source_set = i.set,
                                   source_start = i.el_start, source_end = i.el_end, source_strand = i.el_strand,
                                   source_distance = i.dist)]
  src[!is.na(source_element), source_chrom := src_chrom]
  src[family == "L1" & !is.na(source_element), source_full_length := (source_end - source_start) >= p_int(params$full_length_l1_bp)]
}
if (!is.null(p_str(params$genes)) && any(!is.na(src$src_chrom))) {
  G <- read_elements(params$genes, "gene")
  g <- G[src[!is.na(src_chrom)], on = .(chrom == src_chrom, start <= src_end, end >= src_start), nomatch = 0,
         .(transduction, gene = name)][, .(overlapping_genes = paste(sort(unique(gene)), collapse = ",")), by = transduction]
  src <- merge(src, g, by = "transduction", all.x = TRUE)
}

## group by source: assigned element, else mapped locus within cluster_window
src[, source_locus := NA_character_]
src[!is.na(source_element), source_locus := sprintf("%s:%d-%d(%s)", source_chrom, source_start, source_end, source_strand)]
cw <- p_int(params$cluster_window)
nl <- src[is.na(source_element) & !is.na(src_chrom)][order(src_chrom, src_start)]
if (nrow(nl)) {
  nl[, brk := c(TRUE, src_chrom[-1] != src_chrom[-.N] | src_start[-1] - cummax(src_end)[-.N] > cw)]
  nl[, grp := cumsum(brk)][, locus := sprintf("%s:%d(no_element)", src_chrom[1], src_start[1]), by = grp]
  src[nl, on = "transduction", source_locus := i.locus]
}
un <- src[mapping %in% c("unmapped", "self_only")]
if (nrow(un) > 1 && p_bool(params$run_ava) && !is.null(which_tool("minimap2"))) {
  ufa <- file.path(tempdir(), "unmapped_transductions.fa")
  write_fasta(un$transduction, tseqs[un$transduction, on = "id", seq], ufa)
  upaf <- file.path(tempdir(), "unmapped_ava.paf")
  run_tool(which_tool("minimap2"), c("-c", "-x", params$minimap2_preset, "--secondary=yes", "-X", "-t", threads, ufa, ufa), stdout = upaf)
  A <- read_paf(upaf)
  if (nrow(A)) {
    A <- A[qname != tname & alen > 0 & matches / alen >= p_num(params$min_identity) &
             (qend - qstart) / pmin(qlen, tlen) >= p_num(params$min_qcov)]
    nodes <- un$transduction
    comp <- setNames(seq_along(nodes), nodes)
    a <- match(A$qname, nodes); b <- match(A$tname, nodes)
    repeat {   # label propagation to connected components
      lab <- pmin(comp[a], comp[b])
      upd <- data.table(n = c(a, b), l = c(lab, lab))[, .(l = min(l)), by = n]
      new <- comp; new[upd$n] <- pmin(new[upd$n], upd$l)
      if (identical(new, comp)) break
      comp <- new
    }
    cs <- table(comp)
    multi <- names(cs)[cs > 1]
    if (length(multi)) {
      ids <- match(as.character(comp), multi)
      src[match(nodes[!is.na(ids)], transduction), source_locus := sprintf("nonref_cluster_%04d", ids[!is.na(ids)])]
    }
  }
}
src[, outcome := mapping]
src[!is.na(source_element), outcome := "assigned_element"]
src[stri_startswith_fixed(source_locus, "nonref_cluster") %in% TRUE, outcome := "shared_unmapped"]
save_tab(src, "transduction_sources")

srcsum <- src[!is.na(source_locus), .(n_offspring = uniqueN(variant_id), kinds = paste(sort(unique(kind)), collapse = ","),
  families = paste(sort(unique(family)), collapse = ","), source_element = source_element[1],
  median_td_len = as.numeric(median(td_len)), offspring_chroms = paste(sort(unique(ins_chrom)), collapse = ",")),
  by = source_locus][order(-n_offspring)]
save_tab(srcsum, "source_summary")
