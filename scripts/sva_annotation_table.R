#!/usr/bin/env Rscript
## Builds the per-SVA annotation table read by notebooks/sva_enrichment_permutation.Rmd:
## one row per polymorphic (non-reference) and fixed (reference) SVA with its class,
## allele count/frequency, structural class (MEI structural catalog), repeat-length
## statistics and the value of every genomic annotation.
##
## Outputs (<outdir>/<prefix>.*):
##   sva_annotation_table.tsv.gz  one row per SVA
##   alu_annotation_table.tsv.gz  one row per polymorphic Alu insertion, the same
##                                columns as the SVA table minus the repeat lengths
##   l1_annotation_table.tsv.gz   the same for polymorphic L1 insertions (each only
##                                when its bed and genotype table resolve)
##   sva_annotations.rds          the annotation objects, the functions that look
##                                them up, genome, excluded/included regions,
##                                SVA-gene links and genes. The notebook's
##                                positional null scores RANDOM positions, so it
##                                needs these, not only the per-SVA values.
##   *.bed / *.bedgraph / *.tsv   the tracks as built
##
## Stages:
##   1. Arguments and helpers
##   2. Annotation library (LIB) -- saved into the bundle for the notebook
##   3. Self-checks of the library
##   4. Inputs (local, URL, gs:// or relative to --gcloud-bucket)
##   5. SVA table: loci, class, subfamily, orientation, repeat lengths; and the
##      polymorphic Alu / L1 insertions, collapsed the same way
##   6. Gene-derived tracks (bedtools)
##   7. Replication timing and fork-directionality tracks
##   8. Annotations
##   9. Allele counts and frequencies
##  10. Structural class from the MEI structural catalog
##  11. Outputs

suppressPackageStartupMessages({
  library(argparser)
  library(data.table)
  library(dplyr)
  library(stringi)
})

## ==== 1. Arguments and helpers ======================================================

## Every input is a local file, an http(s) URL, a gs:// URL or a path relative
## to --gcloud-bucket; a bare file name is looked up in --ref-folder of that
## bucket. The defaults are the files notebooks/sva_enrichment_permutation.Rmd
## used before this script existed.
ARGS <- list(
  ## where bucket-relative inputs live
  gcloud_bucket = c("", "Workspace bucket (with or without gs://) that bucket-relative inputs are read from"),
  ref_folder = c("ref", "Bucket folder holding the reference files named by bare file name"),
  ## SVAs and their genotypes
  mei_bed = c("scratch/yulia/l1meaid_intactmei_annotations.bcftoolsmerge.wholegenome.phaseII.FILTERED.ori_intact_anno.SVA.clustered.bed",
              "Clustered polymorphic SVA BED (chrom, start, end, REF, ALT, ID, MEI_type, MEI_subfamily, length, orientation, intact_status, cluster_id)"),
  mei_lengths = c("submissions/14807ece-5021-458a-a5f4-1ab1d7be603f/ProcessMEISVAData/68893a11-d287-4253-a534-4440ddc778cf/call-ProcessMEISVADataTask/SVA_MEI_lengths_tidy.txt.gz",
                  "Tidy polymorphic SVA repeat lengths (sva_id, indiv, length_hexamer, length_VNTR)"),
  ref_table = c("SVA_regions_van_bree.bed", "Reference SVA table (SVA_regions_van_bree.bed or _flanked.txt)"),
  ref_lengths = c("submissions/3346a6a4-d765-4e60-b787-f455f654c33f/ProcessRefSVAData/ccff4f75-4917-45ef-9e74-c421380b6ac1/call-ProcessRefSVADataTask/ref_SVA_length_filtered.txt.gz",
                  "Tidy reference SVA repeat lengths"),
  genotype_table = c("submissions/58e986ab-5cbb-494b-ae6e-d7331d0eb4f4/FilterMEIVCF/af69906c-19c8-47da-a5fd-e2500bbdd15d/call-GzipGTClustered/SVA_gt_tidy_clustered.txt.gz",
                     "Tidy clustered genotypes (CHROM, cluster_id, indiv, hap1, hap2, ...)"),
  sample_metadata = c("gs://prod-drc-broad/longreads/demo_group/lrWGS.phase2.sample.final_set.2025-05-27.tsv",
                      "Sample metadata (participant, sex) for chrX/chrY allele numbers and sex filtering"),
  catalog_per_variant = c("", "MEIStructuralCatalog <prefix>.structure_per_variant.tsv: structural class per SVA"),
  ## polymorphic Alu / L1, same layouts as --mei-bed / --genotype-table
  ## (the files notebooks/mei_sva_internal_var.Rmd and sva_methyl_analysis.Rmd read)
  alu_bed = c("scratch/yulia/l1meaid_intactmei_annotations.bcftoolsmerge.wholegenome.phaseII.FILTERED.ori_intact_anno.ALU.clustered.bed",
              "Clustered polymorphic Alu BED (same 12 columns as --mei-bed)"),
  alu_genotype_table = c("submissions/58e986ab-5cbb-494b-ae6e-d7331d0eb4f4/FilterMEIVCF/1a2a317e-02f4-4675-8bdd-3d88de16fdc4/call-GzipGTClustered/ALU_gt_tidy_clustered.txt.gz",
                         "FilterMEIVCF's tidy clustered Alu genotypes"),
  l1_bed = c("scratch/yulia/l1meaid_intactmei_annotations.bcftoolsmerge.wholegenome.phaseII.FILTERED.ori_intact_anno.LINE.clustered.bed",
             "Clustered polymorphic L1 BED (same 12 columns as --mei-bed)"),
  l1_genotype_table = c("submissions/58e986ab-5cbb-494b-ae6e-d7331d0eb4f4/FilterMEIVCF/ab2ab051-f9b9-4fba-854e-ef5aa5231fa0/call-GzipGTClustered/LINE_gt_tidy_clustered.txt.gz",
                        "FilterMEIVCF's tidy clustered L1 genotypes"),
  ## reference genome and gene annotation (public copies used if missing from --ref-folder)
  genome = c("hg38.chrom.sizes", "Chromosome sizes"),
  reference = c("hg38.fa", "Reference FASTA for GC (or give --gc-bedgraph)"),
  gc_bedgraph = c("", "GC bedGraph instead of --reference"),
  gencode_gtf = c("gencode.v48.basic.annotation.gtf.gz", "GENCODE GTF for gene context, deserts, intron orientation, links"),
  loeuf_table = c("gnomad.v4.1.1.constraint_metrics.tsv.bgz", "gnomAD constraint table (gene_constraint annotation)"),
  ## regions the notebook's null may not sample
  exclude_gap = c("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz", "Assembly gaps"),
  exclude_centromeres = c("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz", "Centromeres"),
  exclude_blacklist = c("https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz", "ENCODE blacklist"),
  exclude_bed = c("", "Extra excluded regions"),
  include_bed = c("", "Restrict the null universe to these regions"),
  ## replication timing and fork directionality
  repliseq_early = c("https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/66a139d7-9bd2-44a0-be29-a26bf4ffe980/4DNFIRKOXCUW.bedGraph.gz", "4DN Repli-seq early fraction"),
  repliseq_late = c("https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/2bdfbe08-50c1-4be2-8b17-f64e2857b429/4DNFI5TMO13R.bedGraph.gz", "4DN Repli-seq late fraction"),
  repliseq_bigwig = c("", "Prebuilt replication-timing track (skips building from early/late)"),
  rfd_bedgraph = c("", "hg38 OK-seq RFD track"),
  okseq_watson = c("", "OK-seq Watson coverage (with --okseq-crick)"),
  okseq_crick = c("", "OK-seq Crick coverage"),
  rfd_hg19_bigwig = c("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5814nnn/GSM5814025/suppl/GSM5814025_Hela_S3_all_cutadapt_sort_delDupl_merge_hmm_RFD_cutoff10_bs1kb_sort.bw", "Published hg19 RFD, lifted to hg38"),
  liftover_chain = c("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19 -> hg38 chain"),
  rfd_bin = c("10000", "RFD bin size when built from OK-seq coverage"),
  rfd_flip = c("false", "Invert the RFD polarity"),
  rfd_min_abs = c("0.1", "|RFD| below this is ambiguous"),
  ucsc_exe_base = c("https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64", "UCSC binaries"),
  ## annotation settings
  protein_coding_only = c("false", "Gene annotations from protein-coding genes only"),
  flanks = c("10000,2000", "Upstream/downstream flank distances"),
  desert_distance = c("500000", "Gene-desert distance"),
  link_window = c("10000", "SVA-gene link window"),
  gene_body_levels = c("CDS,UTR,exon,intron", "Gene-body levels in priority order"),
  gene_background = c("intergenic", "gene_context background level"),
  loeuf_high_max = c("0.6", "LOEUF below this is high constraint"),
  numeric_annotations = c("", "Extra numeric tracks, name=path,..."),
  category_annotations = c("", "Extra categorical BEDs, name=path,..."),
  category_levels = c("", "Priority order per categorical annotation, name=a|b|c,..."),
  category_backgrounds = c("", "Background level per categorical annotation, name=bg,..."),
  gc_window = c("10000", "GC flank window (half on each side)"),
  gc_bin_size = c("500", "GC bin size"),
  point_mode = c("midpoint", "Element anchor: midpoint or start"),
  ## per-SVA statistics
  subfamily_majority = c("0.5", "Majority needed for a cluster's subfamily"),
  min_obs = c("1", "Blank repeat-length statistics of elements seen fewer times"),
  hexamer_bin = c("6", "bp bin for distinct hexamer lengths"),
  vntr_bin = c("24", "bp bin for distinct VNTR lengths"),
  polymorphic_label = c("polymorphic", "Class label of non-reference SVAs"),
  fixed_label = c("fixed", "Class label of reference SVAs"),
  ## output
  cache_dir = c("downloads", "Where URL inputs are downloaded"),
  outdir = c(".", "Output directory"),
  prefix = c("sva", "Output file prefix")
)
ap <- arg_parser("Build the per-SVA annotation table for the SVA enrichment permutation notebook")
for (n in names(ARGS)) ap <- add_argument(ap, paste0("--", gsub("_", "-", n)), help = ARGS[[n]][2], default = ARGS[[n]][1])
argv <- parse_args(ap)
params <- setNames(lapply(names(ARGS), function(n) {
  v <- argv[[n]]
  if (is.null(v) || (length(v) == 1 && is.na(v))) "" else as.character(v)
}), names(ARGS))

p_num <- function(x) as.numeric(x)
p_int <- function(x) as.integer(x)
p_bool <- function(x) isTRUE(as.logical(x)) || identical(tolower(as.character(x)), "true")
p_str <- function(x) if (is.null(x) || all(is.na(x)) || identical(x, "")) NULL else as.character(x)
## "a=1,b=2" -> list(a = "1", b = "2")
p_map <- function(x) {
  x <- p_str(x)
  if (is.null(x)) return(list())
  kv <- stri_split_fixed(stri_split_fixed(x, ",")[[1]], "=", n = 2)
  setNames(lapply(kv, function(z) stri_trim_both(z[2])), vapply(kv, function(z) stri_trim_both(z[1]), ""))
}
p_vec <- function(x, sep = ",") {
  x <- p_str(x)
  if (is.null(x)) return(character())
  stri_trim_both(stri_split_fixed(x, sep)[[1]])
}
PL <- params$polymorphic_label
FX <- params$fixed_label

dir.create(params$outdir, showWarnings = FALSE, recursive = TRUE)
out_path <- function(name, ext = "tsv") file.path(params$outdir, paste0(params$prefix, ".", name, ".", ext))
## unquoted: bedtools reads several of these, and no field contains a tab
write_tsv <- function(d, path, col.names = TRUE) fwrite(d, path, sep = "\t", na = "NA", quote = FALSE, col.names = col.names)
show_tab <- function(x, caption = NULL) {
  if (!is.null(caption)) cat(caption, "\n")
  print(as.data.frame(x))
  invisible(NULL)
}
LOG_T0 <- Sys.time()
LOG <- function(fmt, ...) {
  cat(sprintf("[%s +%6.0fs] %s\n", format(Sys.time(), "%H:%M:%S"),
              as.numeric(difftime(Sys.time(), LOG_T0, units = "secs")), sprintf(fmt, ...)), file = stderr())
}

which_tool <- function(name) if (nzchar(Sys.which(name))) unname(Sys.which(name)) else NULL
run_tool <- function(cmd, args, stdout = "") {
  status <- system2(cmd, args, stdout = stdout, stderr = FALSE)
  if (!identical(as.integer(status), 0L)) stop(sprintf("%s failed (exit %s)", cmd, status))
  invisible(status)
}
run_shell <- function(cmd, out) {
  if (system(paste(cmd, ">", shQuote(out))) != 0) stop("failed: ", cmd)
}

CACHE <- params$cache_dir
WORK <- file.path(tempdir(), "sva_annot_work")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

http_download <- function(url, dest = CACHE) {
  outf <- file.path(dest, basename(sub("[?].*$", "", url)))
  if (file.exists(outf) && file.size(outf) > 0) return(outf)
  message("downloading ", url)
  ok <- tryCatch({
    utils::download.file(url, outf, mode = "wb", quiet = TRUE)
    file.exists(outf) && file.size(outf) > 0
  }, error = function(e) { message("  failed: ", conditionMessage(e)); FALSE })
  if (!ok) { unlink(outf); return(NULL) }
  outf
}

## gs:// URL, or a path relative to --gcloud-bucket (as notebooks/utils.R's
## gs_download does), copied into the cache with `gcloud storage cp`.
gs_resolve <- function(path) {
  if (startsWith(path, "gs://")) return(path)
  bucket <- sub("^gs://", "", sub("/+$", "", params$gcloud_bucket))
  if (!nzchar(bucket)) stop("'", path, "' is relative to the workspace bucket, but --gcloud-bucket is not set")
  paste0("gs://", bucket, "/", sub("^/+", "", path))
}
gs_download <- function(path) {
  url <- gs_resolve(path)
  outf <- file.path(CACHE, basename(url))
  if (file.exists(outf) && file.size(outf) > 0) return(outf)
  message("copying ", url)
  st <- system2("gcloud", c("storage", "cp", shQuote(url), shQuote(outf)))
  if (st != 0 || !file.exists(outf)) stop("gcloud storage cp failed for ", url)
  outf
}

## Resolve one input: a local file is used as is (the WDL localizes File
## inputs), an http(s)/ftp URL is downloaded, anything else is a gs:// URL or
## bucket-relative path. NULL when absent or unreachable, so a missing optional
## track costs only that annotation.
fetch <- function(spec, what = "input", required = FALSE) {
  spec <- p_str(spec)
  if (is.null(spec)) {
    if (required) stop(what, " is required")
    return(NULL)
  }
  out <- tryCatch({
    if (file.exists(spec)) spec
    else if (grepl("^(https?|ftp)://", spec)) http_download(spec)
    else gs_download(spec)
  }, error = function(e) { message(what, ": ", conditionMessage(e)); NULL })
  if (is.null(out) && required) stop("could not read required input ", what, ": ", spec)
  if (is.null(out)) message("skipping ", what, ": ", spec, " could not be read")
  out
}

## A bare file name (no slash, no scheme) lives in --ref-folder of the bucket;
## when it is missing there, the public copy (if any) is downloaded instead.
ref_path <- function(x) {
  x <- p_str(x)
  if (is.null(x) || grepl("[/:]", x) || file.exists(x)) return(x)
  paste0(sub("/+$", "", params$ref_folder), "/", x)
}
fetch_ref <- function(x, what, fallback = NULL, required = FALSE) {
  out <- fetch(ref_path(x), what)
  if (is.null(out) && !is.null(fallback)) {
    message("using the public copy of ", what)
    out <- fetch(fallback, what)
  }
  if (is.null(out) && required) stop("could not obtain required input: ", what)
  out
}
UCSC_CHROM_SIZES <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
UCSC_HG38_FA <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"
GENCODE_GTF_URL <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz"

## UCSC binaries (bigWigToBedGraph, liftOver) are fetched on demand; if one
## cannot be obtained, the track that needed it is skipped.
ucsc_tool <- function(name) {
  local <- which_tool(name)
  if (!is.null(local)) return(local)
  dest <- file.path(CACHE, name)
  if (file.exists(dest) && file.access(dest, 1) == 0) return(dest)
  got <- http_download(paste0(sub("/+$", "", params$ucsc_exe_base), "/", name))
  if (!is.null(got)) Sys.chmod(got, "0755")
  got
}

## Plain bedGraph path for a track that may be bigWig or bedGraph.
as_bedgraph <- function(path, what = "track") {
  if (is.null(path) || !grepl("\\.(bw|bigWig|bigwig)$", path, ignore.case = TRUE)) return(path)
  out <- file.path(CACHE, paste0(basename(path), ".bedgraph"))
  if (file.exists(out) && file.size(out) > 0) return(out)
  exe <- ucsc_tool("bigWigToBedGraph")
  if (is.null(exe)) { message("cannot convert ", what, ": bigWigToBedGraph unavailable; skipping"); return(NULL) }
  ok <- tryCatch({ run_tool(exe, c(path, out)); TRUE },
                 error = function(e) { message("bigWigToBedGraph failed on ", what, ": ", conditionMessage(e)); FALSE })
  if (!ok || !file.exists(out) || file.size(out) == 0) return(NULL)
  out
}

## hg19 -> hg38 for tracks that only exist on the older assembly.
liftover_bedgraph <- function(path, chain, what = "track") {
  if (is.null(path) || is.null(chain)) return(NULL)
  out <- file.path(CACHE, paste0(basename(path), ".hg38.bedgraph"))
  if (file.exists(out) && file.size(out) > 0) return(out)
  exe <- ucsc_tool("liftOver")
  if (is.null(exe)) { message("cannot lift ", what, ": liftOver unavailable"); return(NULL) }
  unmapped <- file.path(CACHE, paste0(basename(path), ".unmapped"))
  raw <- file.path(CACHE, paste0(basename(path), ".lifted.raw"))
  ok <- tryCatch({ run_tool(exe, c("-bedPlus=4", path, chain, raw, unmapped)); TRUE },
                 error = function(e) { message("liftOver failed on ", what, ": ", conditionMessage(e)); FALSE })
  if (!ok || !file.exists(raw) || file.size(raw) == 0) return(NULL)
  run_shell(sprintf("sort -k1,1 -k2,2n %s", shQuote(raw)), out)
  n_un <- if (file.exists(unmapped)) sum(!startsWith(readLines(unmapped, warn = FALSE), "#")) else NA
  message(sprintf("lifted %s to hg38 (%s intervals failed to map)", what, format(n_un, big.mark = ",")))
  out
}

## chrom, start, end, value from a bedGraph (any extra columns dropped)
read_bedgraph <- function(path) {
  d <- fread_any(path, header = FALSE, showProgress = FALSE, select = 1:4)
  setnames(d, c("chrom", "start", "end", "value"))
  d[, chrom := as.character(chrom)][]
}

## ==== 2. Annotation library =========================================================
## Interval algebra, genome coordinates, the annotation constructors and their
## lookups. Kept in one environment, LIB, which is saved into the bundle: the
## notebook loads these exact functions, so its null draws use the lookups that
## built this table. They use data.table and stringi unqualified, and the
## notebook attaches neither, so LIB's parent imports both: the functions then
## resolve those calls wherever the bundle is loaded.
lib_imports <- new.env(parent = globalenv())
for (pkg in c("data.table", "stringi")) {
  for (f in getNamespaceExports(pkg)) assign(f, getExportedValue(pkg, f), envir = lib_imports)
}
LIB <- local({
  fread_any <- function(path, ...) {
    if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), ...) else fread(path, ...)
  }

  ## ---- Interval algebra -------------------------------------------------------

  merge_intervals <- function(d) {
    d <- d[end > start, .(chrom, start, end)]
    if (!nrow(d)) return(d)
    setorder(d, chrom, start, end)
    d[, grp := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = chrom]
    out <- d[, .(start = min(start), end = max(end)), by = .(chrom, grp)]
    out[, grp := NULL][]
  }

  ## a minus b (both become merged first)
  subtract_intervals <- function(a, b) {
    a <- merge_intervals(a)
    if (is.null(b) || !nrow(b)) return(a)
    b <- merge_intervals(b)
    a[, aid := .I]
    setnames(b, c("start", "end"), c("bstart", "bend"))
    setkey(b, chrom, bstart, bend)
    ov <- foverlaps(a, b, by.x = c("chrom", "start", "end"), by.y = c("chrom", "bstart", "bend"),
                    type = "any", nomatch = NA)
    keep <- ov[is.na(bstart), .(chrom, start, end)]
    hit <- ov[!is.na(bstart)]
    if (nrow(hit)) {
      setorder(hit, aid, bstart)
      ## a = [S,E] with disjoint sorted clipped b pieces [s1,e1]..[sk,ek]
      ## -> remainder = [S,s1], [e1,s2], ..., [ek,E]
      gaps <- hit[, {
        s <- pmax(bstart, start[1]); e <- pmin(bend, end[1])
        .(chrom = chrom[1], start = c(start[1], e), end = c(s, end[1]))
      }, by = aid][, .(chrom, start, end)]
      keep <- rbind(keep, gaps)
    }
    setnames(b, c("bstart", "bend"), c("start", "end"))
    merge_intervals(keep[end > start])
  }

  intersect_intervals <- function(a, b) {
    a <- merge_intervals(a); b <- merge_intervals(b)
    if (!nrow(a) || !nrow(b)) return(a[0])
    setnames(b, c("start", "end"), c("bstart", "bend"))
    setkey(b, chrom, bstart, bend)
    ov <- foverlaps(a, b, by.x = c("chrom", "start", "end"), by.y = c("chrom", "bstart", "bend"),
                    type = "any", nomatch = 0)
    out <- ov[, .(chrom, start = pmax(start, bstart), end = pmin(end, bend))]
    setnames(b, c("bstart", "bend"), c("start", "end"))
    merge_intervals(out)
  }

  ## ---- Genome and global coordinates ------------------------------------------

  read_genome <- function(path, keep_chroms = NULL) {
    g <- fread_any(path, header = FALSE, colClasses = list(character = 1))
    setnames(g, 1:2, c("chrom", "size"))
    g <- g[, .(chrom = as.character(chrom), size = as.numeric(size))][size > 0]
    if (!is.null(keep_chroms)) g <- g[chrom %in% keep_chroms]
    setorder(g, chrom)
    g[, offset := cumsum(c(0, head(size, -1)))]
    g[, `:=`(gstart = offset, gend = offset + size)]
    g[]
  }

  gpos <- function(genome, chrom, pos) {
    off <- genome$offset[match(chrom, genome$chrom)]
    off + pos
  }

  ## BED -> global intervals, dropping anything off the genome
  to_global <- function(d, genome) {
    d <- copy(d[chrom %in% genome$chrom])
    off <- genome$offset[match(d$chrom, genome$chrom)]
    sz <- genome$size[match(d$chrom, genome$chrom)]
    d[, `:=`(gs = off + pmax(0, start), ge = off + pmin(sz, end))]
    d[ge > gs]
  }

  ## UCSC table dumps (gap.txt.gz, centromeres.txt.gz) lead with a `bin` column,
  ## so the interval is in fields 2-4 rather than 1-3. Detect and shift rather
  ## than silently reading the bin number as a chromosome.
  read_region_file <- function(path) {
    d <- fread_any(path, header = FALSE, showProgress = FALSE)
    if (ncol(d) >= 4 && is.numeric(d[[1]]) && is.character(d[[2]]) &&
        any(grepl("^(chr|[0-9XYM])", d[[2]][seq_len(min(5, nrow(d)))]))) {
      out <- data.table(chrom = as.character(d[[2]]), start = as.integer(d[[3]]), end = as.integer(d[[4]]))
    } else {
      out <- data.table(chrom = as.character(d[[1]]), start = as.integer(d[[2]]), end = as.integer(d[[3]]))
    }
    out[!is.na(start) & !is.na(end) & end > start]
  }

  read_bed <- function(path, name_col = 4L) {
    d <- fread_any(path, header = FALSE, fill = TRUE, showProgress = FALSE)
    out <- data.table(chrom = as.character(d[[1]]), start = as.integer(d[[2]]), end = as.integer(d[[3]]))
    if (ncol(d) >= name_col) out[, name := as.character(d[[name_col]])]
    out[!is.na(start) & !is.na(end)]
  }

  ## ---- Gene constraint (LOEUF) --------------------------------------------------

  ## The gnomAD v4.1 constraint table, read as mei_sva_internal_var.Rmd does:
  ## MANE Select transcripts of Ensembl genes, LOEUF = lof.oe_ci.upper.
  read_loeuf <- function(path) {
    d <- if (grepl("\\.b?gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), showProgress = FALSE)
         else fread(path, showProgress = FALSE)
    need <- c("gene_id", "mane_select", "lof.oe_ci.upper")
    if (!all(need %in% names(d))) stop(path, ": needs columns ", paste(need, collapse = ", "))
    d <- d[tolower(as.character(mane_select)) == "true" & startsWith(gene_id, "ENSG")]
    d[, .(LOEUF = suppressWarnings(min(as.numeric(lof.oe_ci.upper), na.rm = TRUE))), by = gene_id][is.finite(LOEUF)]
  }

  CONSTRAINT_CLASSES <- c("high_constraint", "low_constraint", "no_LOEUF")

  ## gene_context split by the constraint of the gene each feature belongs to.
  ## `level_files` maps each gene_context level, in its priority order, to the
  ## per-gene BED it is built from and that BED's gene-id column. Every level is
  ## split into high / low / no_LOEUF and the pieces are listed level-major, so
  ## priority resolution picks the gene_context level first and the constraint
  ## within it: the three pieces of a level always add up to that level's
  ## gene_context count. A feature shared by genes of different constraint
  ## resolves to the more constrained one.
  build_constraint_bed <- function(level_files, loeuf, high_max) {
    cls_of <- function(ids) {
      l <- loeuf$LOEUF[match(sub("\\..*$", "", ids), loeuf$gene_id)]
      fifelse(is.na(l), "no_LOEUF", fifelse(l < high_max, "high_constraint", "low_constraint"))
    }
    pieces <- lapply(names(level_files), function(lv) {
      spec <- level_files[[lv]]
      if (!file.exists(spec$path) || file.size(spec$path) == 0) return(NULL)
      d <- fread(spec$path, header = FALSE, select = c(1:3, spec$id_col), showProgress = FALSE,
                 col.names = c("chrom", "start", "end", "gene_id"))
      d[, .(chrom = as.character(chrom), start, end, name = paste0(lv, "_", cls_of(gene_id)))]
    })
    bed <- rbindlist(pieces)
    levels <- as.vector(outer(CONSTRAINT_CLASSES, names(level_files), function(c, l) paste0(l, "_", c)))
    list(bed = bed, levels = levels)
  }

  ## gene_constraint must be an exact subdivision of gene_context.
  check_constraint_split <- function(gc_ann, con_ann, x, ci, half) {
    a <- as.character(lookup_annotation(gc_ann, x, ci, half))
    b <- sub(paste0("_(", paste(CONSTRAINT_CLASSES, collapse = "|"), ")$"), "", as.character(lookup_annotation(con_ann, x, ci, half)))
    list(n = length(a), mismatches = sum(a != b, na.rm = TRUE) + sum(is.na(a) != is.na(b)))
  }

  ## ---- Annotations ------------------------------------------------------------
  ## Every annotation is a list(type, name, ...) with a vectorised lookup over
  ## global positions. Categorical annotations are resolved by priority: the
  ## first level whose intervals contain the point wins, else `background`.

  make_category_annotation <- function(name, bed, genome, levels, background = "intergenic") {
    iv <- lapply(levels, function(lv) {
      g <- to_global(merge_intervals(bed[name == lv]), genome)
      setorder(g, gs)
      list(level = lv, gs = g$gs, ge = g$ge)
    })
    structure(list(type = "category", name = name, levels = c(levels, background),
                   background = background, iv = iv), class = "sva_ann")
  }

  ## Integer genomic span of an element anchored at x with half-width `half`.
  ## A point-like insertion (half = 0) still occupies one base, so every element
  ## has a span of at least 1 bp and the overlap test below behaves sensibly for
  ## both conventions.
  elem_span <- function(x, half) {
    if (is.null(half)) half <- rep(0, length(x))
    s0 <- x - as.integer(floor(half))
    list(s = s0, e = s0 + pmax(1L, as.integer(round(2 * half))))
  }

  ## Does [s, e) overlap any of the sorted, disjoint intervals [gs, ge)?
  ## Only the last interval starting before e can overlap: the intervals are
  ## merged and sorted, so ge is increasing and everything earlier ends earlier,
  ## while everything later starts at or after e.
  overlaps_any <- function(gs, ge, s, e) {
    if (!length(gs)) return(rep(FALSE, length(s)))
    k <- findInterval(e - 1, gs)
    ok <- k > 0
    res <- logical(length(s))
    res[ok] <- ge[k[ok]] > s[ok]
    res
  }

  inside <- function(gs, ge, x) {
    if (!length(gs)) return(rep(FALSE, length(x)))
    i <- findInterval(x, gs)
    ok <- i > 0
    res <- logical(length(x))
    res[ok] <- x[ok] < ge[i[ok]]
    res
  }

  make_numeric_annotation <- function(name, bg, genome) {
    g <- to_global(bg, genome)
    setorder(g, gs)
    structure(list(type = "numeric", name = name, gs = g$gs, ge = g$ge, value = g$value),
              class = "sva_ann")
  }

  ## Leading/lagging-strand correspondence, from an OK-seq RFD track plus each
  ## element's own orientation.
  ##
  ## RFD = (R - L)/(R + L) with R = rightward-moving forks. For a rightward fork
  ## the nascent leading strand runs 5'->3' with the plus strand's sequence, so
  ## its template is the minus strand; the nascent lagging strand carries the
  ## minus strand's sequence, so the LAGGING-STRAND TEMPLATE is the plus strand.
  ## Hence RFD > 0 -> lagging template = "+", RFD < 0 -> lagging template = "-".
  ##
  ## `flip = TRUE` inverts that mapping; it is for a genuinely inverted RFD
  ## track, which the polarity section diagnoses, and not for making the result
  ## come out the expected way round.
  ##
  ## The categories are named for where the ELEMENT'S SENSE STRAND sits, which
  ## is one complement away from how the literature states the bias. TPRT nicks
  ## one strand and extends it with cDNA complementary to the element mRNA, so
  ## the element's sense strand ends up opposite the nicked strand. The
  ## published preference is for nicking the lagging-strand template (Flasch
  ## 2019; Sultana 2019), which therefore shows up HERE as "leading". An
  ## enrichment in `leading` is the expected result; see the interpretation
  ## notes at the end of the report.
  make_strand_match_annotation <- function(name, bg, genome, flip = FALSE, min_abs = 0) {
    num <- make_numeric_annotation(name, bg, genome)
    structure(list(type = "strand_match", name = name, num = num, flip = flip, min_abs = min_abs,
                   levels = c("lagging", "leading", "ambiguous")), class = "sva_ann")
  }

  ## Intron orientation: is the element inserted in the same direction as the
  ## gene whose intron it lies in? This cannot be a plain interval annotation,
  ## because the answer depends on the ELEMENT's orientation as well as on the
  ## locus -- the same intron is "sense" for one insertion and "antisense" for
  ## another. So it is built like strand_match: an interval lookup returning the
  ## GENE's strand, then a comparison against the element's own.
  ##
  ## Regions where a plus-strand gene's intron overlaps a minus-strand gene's
  ## are removed when the BED is built, not resolved by priority. Priority would
  ## silently call every such element "+", which is a directional bias in
  ## exactly the quantity being measured.
  ## `gate` makes this a true SUBDIVISION of one gene_context level rather than
  ## an independent "does it touch an intron" test. Without it the two disagree
  ## for wide elements: a 2 kb reference SVA can overlap both an exon and an
  ## intron interval, which gene_context resolves to `exon` by priority while a
  ## bare overlap test still counts it as intronic. Gating on the resolved label
  ## makes sense + antisense add up to gene_context=intron exactly.
  make_strand_category_annotation <- function(name, bed, genome, background = "not_intronic",
                                              gate = NULL) {
    cat_ann <- make_category_annotation(name, bed, genome, c("+", "-"), background = "none")
    structure(list(type = "strand_category", name = name, cat = cat_ann, background = background,
                   gate = gate, levels = c("sense", "antisense", background)), class = "sva_ann")
  }

  orientation_match <- function(ann, x, chrom_idx, strand, half = NULL) {
    g <- as.character(lookup_annotation(ann$cat, x, chrom_idx, half))
    out <- rep(ann$background, length(x))
    ok <- !is.na(g) & g %in% c("+", "-") & !is.na(strand)
    if (!is.null(ann$gate)) {
      gl <- as.character(lookup_annotation(ann$gate$ann, x, chrom_idx, half))
      ok <- ok & !is.na(gl) & gl == ann$gate$level
    }
    out[ok] <- ifelse(strand[ok] == g[ok], "sense", "antisense")
    factor(out, levels = ann$levels)
  }

  lagging_template <- function(ann, x, chrom_idx) {
    v <- lookup_annotation(ann$num, x, chrom_idx)
    pos_is_plus <- !ann$flip
    out <- rep(NA_character_, length(v))
    ok <- !is.na(v) & abs(v) >= ann$min_abs
    out[ok & v > 0] <- if (pos_is_plus) "+" else "-"
    out[ok & v < 0] <- if (pos_is_plus) "-" else "+"
    out
  }

  ## GC of the sequence FLANKING an element, from per-bin G+C and A+C+G+T counts
  ## (cumulative sums make any window O(1)). Bins are clamped to the chromosome
  ## so a flank near a chromosome end never reads into its neighbour.
  ##
  ## Why flanks and not a window centred on the element: a reference SVA occupies
  ## real reference sequence, and that sequence is strongly GC-rich (the CCCTCT
  ## hexamer, the GC-rich VNTR, SINE-R). A centred window would therefore measure
  ## mostly the element's own composition and manufacture a "fixed SVAs sit in
  ## GC-rich DNA" result out of nothing. A non-reference insertion has no such
  ## sequence in the reference at all, so the two classes would not even be
  ## measuring the same thing. Taking `window/2` of flank on each side of the
  ## element span makes both classes comparable: every element contributes the
  ## same amount of flanking sequence, and none of its own.
  make_gc_annotation <- function(name, gc_bins, genome, window, bin_size) {
    nb <- ceiling(max(genome$gend) / bin_size)
    num <- numeric(nb); den <- numeric(nb)
    i <- (gc_bins$gs %/% bin_size) + 1L
    ok <- i >= 1 & i <= nb
    num[i[ok]] <- gc_bins$num[ok]; den[i[ok]] <- gc_bins$den[ok]
    structure(list(type = "numeric", name = name, gc = TRUE, bin_size = bin_size, window = window,
                   cnum = c(0, cumsum(num)), cden = c(0, cumsum(den)),
                   chrom_lo = genome$gstart %/% bin_size, chrom_hi = (genome$gend - 1) %/% bin_size,
                   genome = genome), class = "sva_ann")
  }

  ## `half` is each element's half-width in bp (0 for an insertion point). It is
  ## only used by GC, which excludes the element's own span.
  lookup_annotation <- function(ann, x, chrom_idx = NULL, half = NULL) {
    if (ann$type == "category") {
      ## ANY OVERLAP of the element's span counts, resolved by priority (first
      ## level in the list wins). This is the same semantics as the notebook's
      ## `bedtools intersect` + `summarize(annot = min(annot))`: a reference SVA
      ## whose edge clips an exon is called exonic, not intronic.
      ##
      ## A wide element overlaps more features than a point does, but that is not
      ## a bias here because the null preserves each element's width -- its null
      ## counterpart overlaps correspondingly more too, so fold enrichment stays
      ## interpretable. (The one place width is NOT absorbed is the
      ## polymorphic-vs-fixed contrast; see the interpretation notes.)
      sp <- elem_span(x, half)
      out <- rep(ann$background, length(x))
      done <- logical(length(x))
      for (lv in ann$iv) {
        todo <- which(!done)
        if (!length(todo)) break
        hit <- overlaps_any(lv$gs, lv$ge, sp$s[todo], sp$e[todo])
        out[todo[hit]] <- lv$level
        done[todo[hit]] <- TRUE
      }
      return(factor(out, levels = ann$levels))
    }
    if (isTRUE(ann$gc)) {
      bs <- ann$bin_size
      nb_flank <- max(1L, as.integer(ann$window %/% 2) %/% bs)   # bins of flank per side
      sp <- elem_span(x, half)
      es <- sp$s           # element start
      ee <- sp$e           # element end
      ## Only whole bins strictly outside the element are used, so no element
      ## sequence can leak into the flank; the cost is up to one bin of flank
      ## per side. Exclusion is therefore at bin resolution -- use a fine GC
      ## track (UCSC gc5Base is 5 bp) when elements are short.
      hi1 <- (es %/% bs) - 1L
      lo1 <- hi1 - nb_flank + 1L
      lo2 <- ((ee - 1) %/% bs) + 1L
      hi2 <- lo2 + nb_flank - 1L
      clo <- ann$chrom_lo[chrom_idx]; chi <- ann$chrom_hi[chrom_idx]
      lo1 <- pmax(lo1, clo); hi1 <- pmin(hi1, chi)
      lo2 <- pmax(lo2, clo); hi2 <- pmin(hi2, chi)
      seg <- function(lo, hi) {
        ok <- hi >= lo
        n <- numeric(length(lo)); d <- numeric(length(lo))
        n[ok] <- ann$cnum[hi[ok] + 2] - ann$cnum[lo[ok] + 1]
        d[ok] <- ann$cden[hi[ok] + 2] - ann$cden[lo[ok] + 1]
        list(n = n, d = d)
      }
      L <- seg(lo1, hi1); R <- seg(lo2, hi2)
      num <- L$n + R$n; den <- L$d + R$d
      return(ifelse(den > 0, num / den, NA_real_))
    }
    i <- findInterval(x, ann$gs)
    v <- rep(NA_real_, length(x))
    ok <- i > 0
    ok[ok] <- x[ok] < ann$ge[i[ok]]
    v[ok] <- ann$value[i[ok]]
    v
  }

  ## Every annotation's value for each element, labelled exactly as the tests
  ## count it: categorical levels with their background, leading / lagging /
  ## ambiguous, sense / antisense / not_intronic, numeric values as measured.
  ## Unlike annotation_values() nothing is blanked to NA for the mark tests.
  element_annotations <- function(anns, x, genome, strand = NULL, half = NULL) {
    ci <- chrom_index(genome, x)
    out <- list()
    for (a in anns) {
      if (a$type == "strand_match") {
        lt <- lagging_template(a, x, ci)
        out[[a$name]] <- as.character(factor(ifelse(is.na(lt) | is.na(strand), "ambiguous",
                                                    ifelse(strand == lt, "lagging", "leading")), levels = a$levels))
        out[[paste0(a$name, "_value")]] <- as.numeric(lookup_annotation(a$num, x, ci))
      } else if (a$type == "strand_category") {
        out[[a$name]] <- as.character(orientation_match(a, x, ci, strand, half))
      } else if (a$type == "category") {
        out[[a$name]] <- as.character(lookup_annotation(a, x, ci, half))
      } else {
        out[[a$name]] <- as.numeric(lookup_annotation(a, x, ci, half))
      }
    }
    as.data.table(out)
  }


  chrom_index <- function(genome, x) findInterval(x, genome$gstart)

  environment()
}, envir = new.env(parent = lib_imports))
list2env(as.list(LIB), envir = environment())

## ==== 3. Self-checks of the library =================================================
## Run on every build before anything is read or written; a failure stops the
## build. The results travel in the bundle and head the notebook's self-check table.
.checks <- list()
ok <- function(lbl, cond) {
  .checks[[length(.checks) + 1]] <<- data.table(check = lbl, result = if (isTRUE(cond)) "PASS" else "FAIL")
  invisible(NULL)
}
dt <- function(...) data.table(...)
local({
## two 1 kb chromosomes
g <- read_genome({ f <- tempfile(fileext = ".txt"); writeLines(c("c1\t1000", "c2\t1000"), f); f })
## merge
m <- merge_intervals(dt(chrom = c("c1","c1","c1","c2"), start = c(0,5,100,0), end = c(10,20,110,5)))
ok("merge: overlapping merged, disjoint kept", nrow(m) == 3 && m[chrom=="c1"][1, end] == 20)
m2 <- merge_intervals(dt(chrom = "c1", start = c(0, 10), end = c(10, 20)))
ok("merge: book-ended intervals merge", nrow(m2) == 1 && m2$end == 20)

## subtract: hole in the middle
s <- subtract_intervals(dt(chrom="c1", start=0, end=100), dt(chrom="c1", start=40, end=60))
ok("subtract: middle hole -> 2 pieces", nrow(s)==2 && all(s$start==c(0,60)) && all(s$end==c(40,100)))
## subtract: b covers a entirely
s2 <- subtract_intervals(dt(chrom="c1", start=10, end=20), dt(chrom="c1", start=0, end=100))
ok("subtract: full cover -> empty", nrow(s2)==0)
## subtract: b outside a
s3 <- subtract_intervals(dt(chrom="c1", start=0, end=10), dt(chrom="c1", start=50, end=60))
ok("subtract: no overlap -> unchanged", nrow(s3)==1 && s3$end==10)
## subtract: different chromosome must not interfere
s4 <- subtract_intervals(dt(chrom="c1", start=0, end=10), dt(chrom="c2", start=0, end=10))
ok("subtract: other chrom ignored", nrow(s4)==1 && s4$chrom=="c1" && s4$end==10)
## subtract: two holes + edge-touching
s5 <- subtract_intervals(dt(chrom="c1", start=0, end=100),
                         dt(chrom="c1", start=c(0,40,90), end=c(10,60,100)))
ok("subtract: edge holes -> interior pieces only",
   nrow(s5)==2 && all(s5$start==c(10,60)) && all(s5$end==c(40,90)))
## total bp conservation
set.seed(3)
a <- merge_intervals(dt(chrom="c1", start=seq(0,9000,by=1000), end=seq(500,9500,by=1000)))
b <- dt(chrom="c1", start=c(250, 3300), end=c(1200, 3400))
s6 <- subtract_intervals(a, b)
lost <- sum(a$end-a$start) - sum(s6$end-s6$start)
expect_lost <- sum(intersect_intervals(a,b)[, end-start])
ok("subtract: bp conservation vs intersect", lost == expect_lost && lost > 0)

## intersect
i1 <- intersect_intervals(dt(chrom="c1", start=0, end=100), dt(chrom="c1", start=50, end=200))
ok("intersect: clipped correctly", nrow(i1)==1 && i1$start==50 && i1$end==100)

## GC window: a pure-GC region must read 1.0, pure-AT 0.0, and clamp at chrom end
gcb <- data.table(gs = c(0:9, 1000:1009) * 100, num = c(rep(100,10), rep(0,10)), den = rep(100, 20))
ga <- make_gc_annotation("gc", gcb, g, window = 300, bin_size = 100)
v <- lookup_annotation(ga, c(500, 1000*100 + 500), chrom_index(g, c(500, 100500)))
ok("gc: GC-rich window reads 1.0", isTRUE(all.equal(v[1], 1)))
vend <- lookup_annotation(ga, 950, chrom_index(g, 950))
ok("gc: window near chrom end stays finite", is.finite(vend))

## ---- leading/lagging strand correspondence ---------------------------------
## RFD > 0 (rightward forks) => lagging-strand TEMPLATE is the plus strand.
rfd <- dt(chrom = c("c1", "c1", "c1"), start = c(0, 300, 600), end = c(300, 600, 900),
          value = c(0.8, -0.8, 0.02))
sm <- make_strand_match_annotation("lag", rfd, g, flip = FALSE, min_abs = 0.1)
xs <- c(100, 400, 700)
lt <- lagging_template(sm, xs, chrom_index(g, xs))
ok("lagging: RFD>0 -> plus is lagging template", identical(lt[1], "+"))
ok("lagging: RFD<0 -> minus is lagging template", identical(lt[2], "-"))
ok("lagging: |RFD| below min_abs -> NA", is.na(lt[3]))
smf <- make_strand_match_annotation("lag", rfd, g, flip = TRUE, min_abs = 0.1)
ok("lagging: flip inverts the mapping",
   identical(lagging_template(smf, xs, chrom_index(g, xs))[1], "-"))

## ---- intron sense / antisense ----------------------------------------------
## The level depends on BOTH the locus and the element's own orientation: the
## same intron is "sense" for one insertion and "antisense" for another.
iob <- dt(chrom = c("c1", "c1"), start = c(0, 400), end = c(200, 600),
          name = c("+", "-"))
ioa <- make_strand_category_annotation("io", iob, g)
xi <- c(100, 100, 500, 500, 800)
si <- c("+", "-", "+", "-", "+")
om <- orientation_match(ioa, xi, chrom_index(g, xi), si)
ok("intron orientation: + element in + intron is sense", identical(as.character(om[1]), "sense"))
ok("intron orientation: - element in + intron is antisense", identical(as.character(om[2]), "antisense"))
ok("intron orientation: + element in - intron is antisense", identical(as.character(om[3]), "antisense"))
ok("intron orientation: - element in - intron is sense", identical(as.character(om[4]), "sense"))
ok("intron orientation: outside any intron is background",
   identical(as.character(om[5]), "not_intronic"))
ok("intron orientation: NA element strand falls back to background",
   identical(as.character(orientation_match(ioa, 100, chrom_index(g, 100), NA_character_)), "not_intronic"))
## gating: with a gene_context gate the split must be a strict subdivision of
## that annotation's `intron` level, not a wider "touches an intron" test
gb <- dt(chrom = c("c1", "c1"), start = c(0, 150), end = c(200, 300),
         name = c("intron", "exon"))
gca <- make_category_annotation("gctx", gb, g, c("exon", "intron"), background = "intergenic")
iog <- make_strand_category_annotation("io", iob, g, gate = list(ann = gca, level = "intron"))
## a WIDE element at 120 spans both the intron and the exon interval; priority
## calls it exonic, so the gated split must NOT count it
wx <- 120; wh <- 60
ok("intron orientation: gate excludes an element priority calls exonic",
   identical(as.character(lookup_annotation(gca, wx, chrom_index(g, wx), wh)), "exon") &&
   identical(as.character(orientation_match(iog, wx, chrom_index(g, wx), "+", wh)), "not_intronic"))
## and a narrow element well inside the intron is still counted
ok("intron orientation: gate keeps a purely intronic element",
   identical(as.character(orientation_match(iog, 50, chrom_index(g, 50), "+", 0)), "sense"))
## ---- flanking GC must ignore the element's own sequence ---------------------
## A reference SVA sits on real reference sequence that is itself GC-rich, so a
## centred window would measure the element, not its neighbourhood. These pin
## the behaviour that prevents that.
## Layout on c1 (bins of 100 bp): bins 0-1 AT, bins 2-3 all-GC, bins 4-5 AT.
gcb2 <- data.table(gs = (0:9) * 100,
                   num = c(0, 0, 100, 100, 0, 0, 0, 0, 0, 0),
                   den = rep(100, 10))
ga2 <- make_gc_annotation("gc", gcb2, g, window = 400, bin_size = 100)
ci2 <- chrom_index(g, 250)
## a point-like insertion at 250 (half = 0) sits inside the GC-rich block and
## its 200 bp flanks per side straddle it
v_point <- lookup_annotation(ga2, 250, ci2, half = 0)
## an element of width 200 centred on 250 covers bins 2-3 exactly; stepping
## over it must leave only AT flanks -> GC 0
v_wide <- lookup_annotation(ga2, 250, ci2, half = 100)
ok("GC: point-like element sees the GC-rich sequence it sits in", v_point > 0)
ok("GC: wide element's own GC-rich span is excluded from its flanks",
   isTRUE(all.equal(v_wide, 0)))
ok("GC: widening an element strictly lowers measured flank GC here", v_wide < v_point)
## with no half supplied the behaviour matches half = 0 (back-compatible)
ok("GC: missing half behaves as a point element",
   isTRUE(all.equal(lookup_annotation(ga2, 250, ci2), v_point)))
## an element wider than the available flank yields NA rather than nonsense
ok("GC: element wider than the chromosome's flank space gives NA",
   is.na(lookup_annotation(ga2, 250, ci2, half = 1e6)))

})
ANNOTATION_CHECKS <- rbindlist(.checks)
n_bad <- sum(ANNOTATION_CHECKS$result != "PASS")
cat(sprintf("Annotation self-checks: %d of %d passed.\n", nrow(ANNOTATION_CHECKS) - n_bad, nrow(ANNOTATION_CHECKS)))
if (n_bad) {
  print(ANNOTATION_CHECKS[result != "PASS"])
  stop("annotation self-check(s) failed; nothing was written")
}

## ==== 4. Inputs =====================================================================
LOG("resolving inputs")
IN <- list(
  mei_bed         = fetch(params$mei_bed, "mei_bed", required = TRUE),
  mei_lengths     = fetch(params$mei_lengths, "mei_lengths", required = TRUE),
  ref_table       = fetch_ref(params$ref_table, "ref_table", required = TRUE),
  ref_lengths     = fetch(params$ref_lengths, "ref_lengths", required = TRUE),
  genotype_table  = fetch(params$genotype_table, "genotype_table", required = TRUE),
  sample_metadata = fetch(params$sample_metadata, "sample_metadata"),
  catalog         = fetch(params$catalog_per_variant, "catalog_per_variant"),
  genome          = fetch_ref(params$genome, "genome", fallback = UCSC_CHROM_SIZES, required = TRUE),
  reference       = fetch_ref(params$reference, "reference", fallback = UCSC_HG38_FA),
  gc_bedgraph     = as_bedgraph(fetch(params$gc_bedgraph, "gc_bedgraph"), "gc_bedgraph"),
  gencode_gtf     = fetch_ref(params$gencode_gtf, "gencode_gtf", fallback = GENCODE_GTF_URL),
  loeuf           = fetch_ref(params$loeuf_table, "loeuf_table"),
  include_bed     = fetch(params$include_bed, "include_bed"),
  rfd             = as_bedgraph(fetch(params$rfd_bedgraph, "rfd_bedgraph"), "rfd")
)
## The other polymorphic MEI families, each with its own element ids (as the
## notebooks name them) and the names its insertions carry in the structural
## catalog's `family` column.
OTHER_FAMILIES <- list(
  Alu = list(bed = fetch(params$alu_bed, "alu_bed"), genotype_table = fetch(params$alu_genotype_table, "alu_genotype_table"),
             id_fmt = "Alu_MEI_cl%05d", catalog_families = "ALU"),
  L1 = list(bed = fetch(params$l1_bed, "l1_bed"), genotype_table = fetch(params$l1_genotype_table, "l1_genotype_table"),
            id_fmt = "LINE1_MEI_cl%05d", catalog_families = c("L1", "LINE1", "LINE"))
)
for (fam in names(OTHER_FAMILIES)) {
  f <- OTHER_FAMILIES[[fam]]
  if (xor(is.null(f$bed), is.null(f$genotype_table))) {
    cat(sprintf("WARNING: %s needs both its bed and its genotype table; it is skipped.\n", fam))
  }
}
OTHER_FAMILIES <- Filter(function(f) !is.null(f$bed) && !is.null(f$genotype_table), OTHER_FAMILIES)
## Extra categorical (BED, level in column 4) and numeric (bedGraph) tracks,
## given as name=path pairs.
IN$category <- Filter(Negate(is.null), lapply(p_map(params$category_annotations), fetch, what = "category annotation"))
IN$numeric <- lapply(p_map(params$numeric_annotations), fetch, what = "numeric annotation")

## hg38 models centromeres as sequence rather than as gaps, so unlike hg19 the
## gap table does not contain them and they are excluded separately.
excl_parts <- Filter(Negate(is.null), list(
  gap         = fetch(params$exclude_gap, "assembly gaps"),
  centromeres = fetch(params$exclude_centromeres, "centromeres"),
  blacklist   = fetch(params$exclude_blacklist, "ENCODE blacklist"),
  extra       = fetch(params$exclude_bed, "exclude_bed")
))
excl <- if (length(excl_parts)) merge_intervals(rbindlist(lapply(excl_parts, read_region_file), fill = TRUE)) else NULL
incl <- if (!is.null(IN$include_bed)) read_bed(IN$include_bed) else NULL

## ==== 5. SVA table ==================================================================
## One row per SVA, reproducing the collapsing in notebooks/mei_sva_internal_var.Rmd
## so loci, orientations and repeat lengths match that notebook.
LOG("building the SVA table")
MEI_COLS <- c("chrom", "start", "end", "REF", "ALT", "ID", "MEI_type", "MEI_subfamily",
              "length", "orientation", "intact_status", "cluster_id")
LENGTH_STATS <- c("mean_hexamer", "max_hexamer", "sd_hexamer", "nuniq_hexamer",
                  "mean_VNTR", "max_VNTR", "sd_VNTR", "nuniq_VNTR")

## Elements lost at each stage, per class: a class that silently reaches zero
## looks exactly like one never supplied, and this table tells them apart.
FUNNEL <- tibble(stage = character(), class = character(), n = integer())
funnel <- function(d, stage) {
  counts <- count(d, class, name = "n") %>% mutate(stage = stage, n = as.integer(n))  # forces d (and any earlier stage) first
  FUNNEL <<- bind_rows(FUNNEL, counts)
  invisible(d)
}

## Distinct lengths after flooring to `bin` bp (one repeat unit), so base-level
## sequencing or alignment jitter does not count as a new allele.
n_distinct_bins <- function(x, bin) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_integer_ else length(unique(floor(x / bin)))
}
max_or_na <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  if (length(x)) max(x) else NA_real_
}
stopifnot(n_distinct_bins(c(100, 101, 105, 106, NA), 6) == 2L,
          is.na(n_distinct_bins(c(NA_real_, NA), 6)),
          n_distinct_bins(c(0, 23, 24, 47, 48), 24) == 3L)

## Per-element repeat-length statistics from a tidy length table (one row per
## haplotype carrying the element). n_hexamer / n_VNTR are the haplotypes with a
## measured length: the sample size of each statistic, and the allele count the
## notebook's mark permutation stratifies on. No allele-frequency cutoff by
## default (min_obs = 1): dropping rare elements would bias the positional null
## toward common loci.
length_stats <- function(path, min_obs) {
  d <- fread_any(path, sep = "\t", showProgress = FALSE)
  miss <- setdiff(c("sva_id", "length_hexamer", "length_VNTR"), names(d))
  if (length(miss)) stop(path, ": missing column(s) ", paste(miss, collapse = ", "))
  if (!"indiv" %in% names(d)) d[, indiv := NA]
  bh <- p_num(params$hexamer_bin); bv <- p_num(params$vntr_bin)
  s <- as_tibble(d) %>%
    group_by(sva_id) %>%
    summarise(n_obs = n(),
              n_indiv = if (all(is.na(indiv))) NA_integer_ else n_distinct(indiv),
              mean_hexamer = mean(length_hexamer, na.rm = TRUE),
              max_hexamer = max_or_na(length_hexamer),
              sd_hexamer = sd(length_hexamer, na.rm = TRUE),
              nuniq_hexamer = n_distinct_bins(length_hexamer, bh),
              n_hexamer = sum(!is.na(length_hexamer)),
              mean_VNTR = mean(length_VNTR, na.rm = TRUE),
              max_VNTR = max_or_na(length_VNTR),
              sd_VNTR = sd(length_VNTR, na.rm = TRUE),
              nuniq_VNTR = n_distinct_bins(length_VNTR, bv),
              n_VNTR = sum(!is.na(length_VNTR)),
              .groups = "drop")
  if (min_obs > 1) s <- s %>% mutate(across(all_of(LENGTH_STATS), ~ if_else(n_obs < min_obs, NA_real_, as.numeric(.x))))
  s
}

## Polymorphic insertions of one family: one locus per cluster (chrom, min
## start, max end), the majority subfamily (clusters without a >
## subfamily_majority majority are dropped, the notebook's ambig_subfamily
## filter) and the majority orientation ("NONE" stays NA; a tie goes to the
## orientation listed first). `lengths` (SVA only) adds the repeat lengths, and
## only the SVA call feeds the element funnel.
read_mei_svas <- function(path, lengths = NULL, id_fmt = "SVA_MEI_cl%05d") {
  b <- fread_any(path, sep = "\t", header = FALSE, showProgress = FALSE)
  if (ncol(b) < length(MEI_COLS)) {
    stop(path, ": expected ", length(MEI_COLS), " columns (", paste(MEI_COLS, collapse = ", "), "), found ", ncol(b))
  }
  setnames(b, seq_along(MEI_COLS), MEI_COLS)
  b <- as_tibble(b) %>%
    mutate(chrom = as.character(chrom), cluster_id = as.integer(cluster_id), row = row_number()) %>%
    filter(!is.na(cluster_id)) %>%
    mutate(sva_id = sprintf(id_fmt, cluster_id))

  loci <- b %>%
    group_by(cluster_id, sva_id) %>%
    summarise(chrom = first(chrom), start = min(start), end = max(end),
              ins_length = as.numeric(median(suppressWarnings(as.numeric(length)), na.rm = TRUE)),
              n_alleles = n(), n_chrom = n_distinct(chrom), .groups = "drop")
  if (any(loci$n_chrom > 1)) cat("Some clusters span more than one chromosome; the first chromosome was kept.\n")

  maj <- p_num(params$subfamily_majority)
  subfamily <- b %>%
    count(cluster_id, MEI_subfamily) %>%
    group_by(cluster_id) %>%
    filter(max(n / sum(n)) > maj) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(cluster_id, subfamily = MEI_subfamily)
  cat(sprintf("%d MEI clusters, %d with an unambiguous subfamily (>%.0f%% majority).\n",
              n_distinct(b$cluster_id), nrow(subfamily), 100 * maj))

  orientation <- b %>%
    filter(orientation %in% c("+", "-")) %>%
    group_by(cluster_id, orientation) %>%
    summarise(n = n(), first_row = min(row), .groups = "drop") %>%
    arrange(cluster_id, desc(n), first_row) %>%
    distinct(cluster_id, .keep_all = TRUE) %>%
    select(cluster_id, orientation)

  svas <- loci %>%
    select(-n_chrom) %>%
    left_join(subfamily, by = "cluster_id") %>%
    left_join(orientation, by = "cluster_id") %>%
    mutate(class = PL)
  if (is.null(lengths)) return(svas)
  funnel(svas, "1. read from its input file")
  funnel(svas, "2. region parsed to coordinates")
  svas <- svas %>% left_join(lengths, by = "sva_id")
  funnel(svas, "3. joined to ref_lengths")
}

## Reference SVAs. Two layouts are in circulation and are NOT interchangeable,
## so the layout is detected rather than assumed:
##   A. SVA_regions_van_bree_flanked.txt -- "chrom:start-end" in column 1, then
##      subfamily, strand, score, SVA_id; spans carry a 50 bp flank each side.
##   B. SVA_regions_van_bree.bed -- chrom, start, end, sva_id, then strand,
##      score, subfamily. Its `#` header says "score, strand" while the data is
##      strand then score; trusting it empties every orientation-based
##      annotation without failing. So columns are identified by CONTENT:
##      strand is the column that is all +/-/., subfamily the one with SVA_* names.
read_ref_svas <- function(path, lengths) {
  r <- fread_any(path, sep = "\t", header = FALSE, showProgress = FALSE)
  if (nrow(r) && grepl("^#", as.character(r[[1]][1]))) r <- r[-1]
  cat(sprintf("ref_table %s: %d rows, %d columns.\n", path, nrow(r), ncol(r)))
  funnel(tibble(class = rep(FX, nrow(r))), "1. read from its input file")
  num <- function(v) suppressWarnings(as.numeric(as.character(v)))

  if (nrow(r) && grepl("^[^:]+:\\d+-\\d+$", as.character(r[[1]][1]))) {
    cat("ref_table layout: region strings ('chrom:start-end' in column 1).\n")
    k <- min(5, ncol(r))
    setnames(r, 1:k, c("region", "subfamily", "strand", "score", "SVA_id")[1:k])
    mm <- stri_match_first_regex(r$region, "^([^:]+):(\\d+)-(\\d+)$")
    rr <- tibble(sva_id = r$SVA_id, chrom = mm[, 2], start = num(mm[, 3]), end = num(mm[, 4]),
                 subfamily = r$subfamily,
                 orientation = if ("strand" %in% names(r)) r$strand else NA_character_)
  } else {
    if (ncol(r) < 4) stop("ref_table looks like a BED but has only ", ncol(r), " columns; need chrom, start, end, id")
    extra <- setdiff(seq_len(ncol(r)), 1:4)
    first_col <- function(test) {
      hit <- extra[vapply(extra, function(i) test(as.character(r[[i]])), logical(1))]
      if (length(hit)) hit[1] else NA_integer_
    }
    i_strand <- first_col(function(v) all(v %in% c("+", "-", ".")))
    i_sub <- first_col(function(v) any(grepl("^SVA[_-]", v)))
    cat(sprintf("ref_table layout: BED; strand = column %s, subfamily = column %s (found by content).\n",
                ifelse(is.na(i_strand), "none found", i_strand), ifelse(is.na(i_sub), "none found", i_sub)))
    if (is.na(i_strand)) {
      cat("WARNING: no strand column in ref_table: fixed SVAs get no orientation, so their",
          "leading/lagging and intron-orientation values will all be missing.\n")
    }
    rr <- tibble(sva_id = as.character(r[[4]]), chrom = as.character(r[[1]]),
                 start = num(r[[2]]), end = num(r[[3]]),
                 subfamily = if (is.na(i_sub)) NA_character_ else as.character(r[[i_sub]]),
                 orientation = if (is.na(i_strand)) NA_character_ else as.character(r[[i_strand]]))
  }
  n_bad <- sum(is.na(rr$start) | is.na(rr$end))
  if (n_bad) cat(sprintf("%d reference rows had unparseable coordinates and were dropped.\n", n_bad))
  rr <- rr %>%
    filter(!is.na(start), !is.na(end)) %>%
    mutate(orientation = if_else(orientation %in% c("+", "-"), orientation, NA_character_),
           class = FX, ins_length = end - start, n_alleles = NA_integer_)
  funnel(rr, "2. region parsed to coordinates")
  cat(sprintf("ref_table chromosome names look like: %s\n", paste(head(sort(unique(rr$chrom)), 4), collapse = ", ")))
  matched <- sum(rr$sva_id %in% lengths$sva_id)
  cat(sprintf("ref_lengths: %d of %d reference IDs matched (%.1f%%).%s\n", matched, nrow(rr),
              100 * matched / max(1, nrow(rr)),
              if (matched == 0) " IDs differ between ref_table and ref_lengths: no fixed SVA will have repeat lengths." else ""))
  rr <- rr %>% left_join(lengths, by = "sva_id")
  funnel(rr, "3. joined to ref_lengths")
}

min_obs <- p_int(params$min_obs)
mei_len <- length_stats(IN$mei_lengths, min_obs)
ref_len <- length_stats(IN$ref_lengths, min_obs)
## A reference element that is also a called MEI is dropped, as in the
## notebook, so the two classes stay disjoint.
n_both <- sum(ref_len$sva_id %in% mei_len$sva_id)
ref_len <- ref_len %>% filter(!sva_id %in% mei_len$sva_id)
if (n_both) cat(sprintf("%d reference elements dropped for also being called MEIs.\n", n_both))

KEEP <- c("chrom", "start", "end", "name", "class", "subfamily", "orientation",
          "n_obs", "n_indiv", "n_alleles", "ins_length", "n_hexamer", "n_VNTR", LENGTH_STATS)
S <- bind_rows(read_mei_svas(IN$mei_bed, mei_len), read_ref_svas(IN$ref_table, ref_len)) %>%
  rename(name = sva_id) %>%
  funnel("4. after dropping elements also called MEIs") %>%
  filter(!is.na(chrom), !is.na(start), !is.na(end)) %>%
  funnel("5. with usable coordinates") %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  arrange(chrom, start) %>%
  select(any_of(KEEP))

## polymorphic Alu / L1, collapsed like the SVAs; `family` names the family
OTHER <- bind_rows(lapply(names(OTHER_FAMILIES), function(fam) {
  cat(sprintf("%s: ", fam))
  read_mei_svas(OTHER_FAMILIES[[fam]]$bed, id_fmt = OTHER_FAMILIES[[fam]]$id_fmt) %>%
    rename(name = sva_id) %>%
    filter(!is.na(chrom), !is.na(start), !is.na(end)) %>%
    mutate(family = fam, start = as.integer(start), end = as.integer(end))
}))

genome <- read_genome(IN$genome, keep_chroms = unique(c(S$chrom, OTHER$chrom)))
## A naming mismatch ("1" vs "chr1") would silently remove a whole class.
lost <- S %>% filter(!chrom %in% genome$chrom) %>% count(class)
for (i in seq_len(nrow(lost))) {
  cat(sprintf("WARNING: %d %s elements are on chromosomes absent from the genome file (which has %s ...).\n",
              lost$n[i], lost$class[i], paste(head(sort(genome$chrom), 4), collapse = ", ")))
}

## The two classes have different interval semantics: a polymorphic SVA is an
## insertion POINT (its sequence is absent from the reference), a fixed SVA
## spans its own 1-3 kb of reference sequence. `point` anchors each element,
## `half` is its half-width, which the GC lookup steps over.
S <- S %>%
  filter(chrom %in% genome$chrom) %>%
  funnel("6. on a chromosome in the genome file") %>%
  mutate(chrom = as.character(chrom), class = as.character(class), subfamily = as.character(subfamily),
         strand = if_else(orientation %in% c("+", "-"), orientation, NA_character_),
         width = pmax(0, as.numeric(end) - as.numeric(start)),
         point = if (params$point_mode == "start") as.numeric(start) else floor((as.numeric(start) + as.numeric(end)) / 2),
         half = width / 2,
         g = gpos(genome, chrom, point))
cat(sprintf("%d SVAs on %d chromosomes.\n", nrow(S), n_distinct(S$chrom)))
show_tab(S %>% group_by(class) %>%
           summarise(n = n(), median_width = as.numeric(median(width)), max_width = max(width),
                     pct_pointlike = round(100 * mean(width <= 1), 1)),
         caption = "Interval conventions by class (point-like = width <= 1 bp)")
if (nrow(OTHER)) {
  OTHER <- OTHER %>%
    filter(chrom %in% genome$chrom) %>%
    mutate(strand = if_else(orientation %in% c("+", "-"), orientation, NA_character_),
           width = pmax(0, as.numeric(end) - as.numeric(start)),
           point = if (params$point_mode == "start") as.numeric(start) else floor((as.numeric(start) + as.numeric(end)) / 2),
           half = width / 2,
           g = gpos(genome, chrom, point)) %>%
    arrange(chrom, start)
  show_tab(count(OTHER, family), caption = "Polymorphic Alu / L1 insertions on the genome's chromosomes")
}

## ==== 6. Gene-derived tracks (bedtools) =============================================
## Gene context, gene deserts, strand-aware introns and SVA-gene links, all from
## the GTF. The shell script writes BEDs into WORK; the R below reads them back
## and checks each, since the script stops at its first failure.
flank_sizes <- as.integer(p_vec(params$flanks))
flank_sizes <- sort(flank_sizes[!is.na(flank_sizes) & flank_sizes > 0])
built_cat <- list()        # name -> list(path, levels, background)
built_strandcat <- list()  # name -> path; the level is the host gene's strand
links <- NULL
genes_named <- NULL

GENE_TRACKS_SH <- r"---(
set -euo pipefail
cd "$AW_WORK"

if [[ "$AW_GTF" == *.gz ]]; then CAT="gzip -dc"; else CAT="cat"; fi

# One GTF pass per feature type, keeping only the analysis chromosomes and
# (optionally) protein-coding genes; column 4 is the gene id.
extract() {                       # $1 = regex matching the feature column
  $CAT "$AW_GTF" | awk -v PC="$AW_PC" -v FEAT="$1" -F'\t' 'BEGIN{OFS="\t"}
      NR==FNR            {ok[$1]=1; next}
      !($1 in ok)        {next}
      PC==1 && $0 !~ /gene_type "protein_coding"/ {next}
      $3 ~ FEAT          {id=$0; sub(/.*gene_id "/, "", id); sub(/".*/, "", id); print $1, $4-1, $5, id}' chroms.txt - \
    | sort -k1,1 -k2,2n
}
extract '^gene$'                                    > genes.bed
extract '^exon$'                                    > exons.bed
extract '^CDS$'                                     > cds.bed
extract '^(UTR|five_prime_UTR|three_prime_UTR)$'    > utr.bed

# Gene bodies with a symbol and a strand, for the flanks and the SVA-gene links
$CAT "$AW_GTF" | awk -v PC="$AW_PC" -F'\t' 'BEGIN{OFS="\t"}
    NR==FNR     {ok[$1]=1; next}
    !($1 in ok) {next}
    PC==1 && $0 !~ /gene_type "protein_coding"/ {next}
    $3=="gene"  {
      id=$0; sub(/.*gene_id "/, "", id); sub(/".*/, "", id)
      nm=$0; if (nm ~ /gene_name "/) { sub(/.*gene_name "/, "", nm); sub(/".*/, "", nm) } else { nm=id }
      print $1, $4-1, $5, nm, id, $7
    }' chroms.txt - | sort -k1,1 -k2,2n > genes6.bed

if [[ ! -s genes6.bed ]]; then
  echo "no genes survived the chromosome filter: the GTF's chromosome naming probably does not match the genome file (GENCODE and UCSC use 'chr1', Ensembl uses '1')" >&2
  exit 1
fi
echo "genes on the analysis chromosomes: $(wc -l < genes6.bed)"

# introns = gene bodies minus exons
"$AW_BEDTOOLS" subtract -a genes.bed -b exons.bed > introns.bed

# Strand-aware introns for the sense/antisense split. Where a + gene's intron
# overlaps a - gene's, the host strand is undefined, so the overlap is cut out
# of both rather than assigned to one.
"$AW_BEDTOOLS" subtract -a genes6.bed -b exons.bed | sort -k1,1 -k2,2n > introns6.bed
awk -F'\t' 'BEGIN{OFS="\t"} $6=="+"' introns6.bed > intron_p.bed
awk -F'\t' 'BEGIN{OFS="\t"} $6=="-"' introns6.bed > intron_m.bed
"$AW_BEDTOOLS" subtract -a intron_p.bed -b intron_m.bed | cut -f1-3 > intron_plus_only.bed
"$AW_BEDTOOLS" subtract -a intron_m.bed -b intron_p.bed | cut -f1-3 > intron_minus_only.bed
echo "strand-unambiguous introns: + $(wc -l < intron_plus_only.bed), - $(wc -l < intron_minus_only.bed)"

# strand-aware flanks, one pair per requested distance
for FL in $AW_FLANKS; do
  "$AW_BEDTOOLS" flank -s -l "$FL" -r 0 -i genes6.bed -g genome.txt > "up_$FL.bed"
  "$AW_BEDTOOLS" flank -s -l 0 -r "$FL" -i genes6.bed -g genome.txt > "down_$FL.bed"
done

# gene deserts: everything further than AW_DESERT from any gene
"$AW_BEDTOOLS" slop -b "$AW_DESERT" -i genes.bed -g genome.txt \
  | "$AW_BEDTOOLS" merge -i - \
  | "$AW_BEDTOOLS" complement -i - -g genome.txt > deserts.bed

# SVA -> gene links within AW_LINKW
"$AW_BEDTOOLS" window -a sva_sorted.bed -b genes6.bed -w "$AW_LINKW" > links.txt
echo "SVA-gene links: $(wc -l < links.txt)"
)---"

work_file <- function(n) file.path(WORK, n)
has_file <- function(n) file.exists(work_file(n)) && file.size(work_file(n)) > 0
## chrom, start, end of a WORK file, labelled `lab`; NULL if absent or empty
read_labelled <- function(n, lab) {
  if (!has_file(n)) return(NULL)
  fread(work_file(n), header = FALSE, select = 1:3, col.names = c("chrom", "start", "end"),
        showProgress = FALSE)[, name := lab][]
}
flank_label <- function(fl) if (fl >= 1000) paste0(fl %/% 1000, "kb") else paste0(fl, "bp")

if (!is.null(IN$gencode_gtf)) {
  bt <- which_tool("bedtools")
  if (is.null(bt)) stop("--gencode-gtf needs bedtools on PATH")
  LOG("building gene-derived tracks")

  ## bedtools flank/slop abort on a chromosome missing from the genome file, and
  ## a GENCODE GTF covers scaffolds, alts and chrM, so the script filters to
  ## the analysis chromosomes. Sizes are written as plain integers: fwrite
  ## would print a round one as "2e+06", which bedtools rejects.
  writeLines(genome$chrom, work_file("chroms.txt"))
  write_tsv(genome[, .(chrom, size = sprintf("%.0f", size))], work_file("genome.txt"), col.names = FALSE)
  write_tsv(bind_rows(S, OTHER) %>% select(chrom, start, end, name) %>% arrange(chrom, start), work_file("sva_sorted.bed"),
            col.names = FALSE)
  Sys.setenv(AW_WORK = WORK,
             AW_GTF = normalizePath(IN$gencode_gtf, mustWork = TRUE),
             AW_BEDTOOLS = bt,
             AW_PC = if (p_bool(params$protein_coding_only)) "1" else "0",
             AW_FLANKS = paste(flank_sizes, collapse = " "),
             AW_DESERT = as.character(p_int(params$desert_distance)),
             AW_LINKW = as.character(p_int(params$link_window)))
  writeLines(GENE_TRACKS_SH, work_file("gene_tracks.sh"))
  status <- system2("bash", work_file("gene_tracks.sh"))
  if (!identical(as.integer(status), 0L)) LOG("WARNING: the gene-track shell step exited %s", status)

  expected <- c("genes.bed", "exons.bed", "cds.bed", "utr.bed", "introns.bed", "genes6.bed",
                "intron_plus_only.bed", "intron_minus_only.bed", "deserts.bed", "links.txt",
                sprintf("up_%d.bed", flank_sizes), sprintf("down_%d.bed", flank_sizes))
  missing_files <- expected[!vapply(expected, has_file, logical(1))]
  if (length(missing_files)) {
    cat(sprintf("WARNING: gene-track step incomplete, missing %s. %s\n", paste(missing_files, collapse = ", "),
                if ("genes6.bed" %in% missing_files) "Without genes6.bed, every gene-based annotation and the SVA-gene links are skipped."
                else "The gene-based annotations are built from whatever did build."))
  }

  ## ONE gene-context annotation. Body features do not depend on the flank
  ## distance, so one annotation per distance would test the same CDS and intron
  ## hypotheses twice. Flanks are listed smallest first and priority resolves
  ## each element to a single label, so the distances come out DISJOINT:
  ## upstream_10kb means "within 10 kb but not within 2 kb".
  flank_levels <- as.vector(rbind(paste0("upstream_", sapply(flank_sizes, flank_label)),
                                  paste0("downstream_", sapply(flank_sizes, flank_label))))
  context_bed <- rbindlist(c(
    list(read_labelled("cds.bed", "CDS"), read_labelled("utr.bed", "UTR"),
         read_labelled("exons.bed", "exon"), read_labelled("introns.bed", "intron")),
    unlist(lapply(flank_sizes, function(fl) list(read_labelled(sprintf("up_%d.bed", fl), paste0("upstream_", flank_label(fl))),
                                                 read_labelled(sprintf("down_%d.bed", fl), paste0("downstream_", flank_label(fl))))),
           recursive = FALSE)), use.names = TRUE)
  if (nrow(context_bed)) {
    setorder(context_bed, chrom, start)
    pth <- out_path("gene_context", "bed")
    write_tsv(context_bed, pth, col.names = FALSE)
    built_cat$gene_context <- list(path = pth, levels = c(p_vec(params$gene_body_levels), flank_levels),
                                   background = params$gene_background)
    cat(sprintf("gene_context: %d intervals (%s)\n", nrow(context_bed),
                paste(sprintf("%s=%d", names(table(context_bed$name)), table(context_bed$name)), collapse = ", ")))
  }

  ## Gene constraint: gene_context's levels, each split by the LOEUF of the gene
  ## the feature belongs to. Its own annotation, so gene_context's counts stay.
  if (nrow(context_bed) && !is.null(IN$loeuf)) {
    level_files <- list(CDS = list(path = work_file("cds.bed"), id_col = 4L),
                        UTR = list(path = work_file("utr.bed"), id_col = 4L),
                        exon = list(path = work_file("exons.bed"), id_col = 4L),
                        intron = list(path = work_file("introns6.bed"), id_col = 5L))
    level_files <- level_files[intersect(p_vec(params$gene_body_levels), names(level_files))]
    for (fl in flank_sizes) {
      level_files[[paste0("upstream_", flank_label(fl))]] <- list(path = work_file(sprintf("up_%d.bed", fl)), id_col = 5L)
      level_files[[paste0("downstream_", flank_label(fl))]] <- list(path = work_file(sprintf("down_%d.bed", fl)), id_col = 5L)
    }
    loeuf <- tryCatch(read_loeuf(IN$loeuf), error = function(e) {
      cat("WARNING: could not read --loeuf-table:", conditionMessage(e), "-- gene_constraint skipped.\n"); NULL })
    if (!is.null(loeuf)) {
      high_max <- p_num(params$loeuf_high_max)
      gcon <- build_constraint_bed(level_files, loeuf, high_max)
      if (nrow(gcon$bed)) {
        setorder(gcon$bed, chrom, start)
        pth <- out_path("gene_constraint", "bed")
        write_tsv(gcon$bed, pth, col.names = FALSE)
        built_cat$gene_constraint <- list(path = pth, levels = gcon$levels, background = params$gene_background)
        gene_ids <- unique(fread(work_file("genes6.bed"), header = FALSE, select = 5L)[[1]])
        gl <- loeuf$LOEUF[match(sub("\\..*$", "", gene_ids), loeuf$gene_id)]
        cl <- table(factor(case_when(is.na(gl) ~ "no_LOEUF", gl < high_max ~ "high_constraint", TRUE ~ "low_constraint"),
                           levels = CONSTRAINT_CLASSES))
        cat(sprintf("gene_constraint: %d MANE genes with LOEUF; annotated genes: %s (high = LOEUF < %s)\n",
                    nrow(loeuf), paste(sprintf("%s %d", names(cl), cl), collapse = ", "), high_max))
      }
    }
  } else if (is.null(IN$loeuf)) {
    cat("No --loeuf-table; the gene_constraint annotation is skipped.\n")
  }

  deserts <- read_labelled("deserts.bed", "gene_desert")
  if (!is.null(deserts) && nrow(deserts)) {
    pth <- out_path("gene_deserts", "bed")
    write_tsv(deserts, pth, col.names = FALSE)
    built_cat$gene_desert <- list(path = pth, levels = "gene_desert", background = "near_gene")
    cat(sprintf("gene_desert: %d deserts (> %s kb from any gene), %.1f%% of the genome\n", nrow(deserts),
                format(p_int(params$desert_distance) / 1000, big.mark = ","),
                100 * sum(deserts$end - deserts$start) / sum(as.numeric(genome$size))))
  } else {
    cat(sprintf("No gene deserts at --desert-distance %s; that annotation is skipped.\n", params$desert_distance))
  }

  ## Intron orientation is its own annotation rather than extra gene_context
  ## levels: gene_context resolves to one label per element, so splitting
  ## `intron` there would lose the intronic count the split is compared against.
  intron_bed <- rbindlist(list(read_labelled("intron_plus_only.bed", "+"),
                               read_labelled("intron_minus_only.bed", "-")), use.names = TRUE)
  if (nrow(intron_bed)) {
    setorder(intron_bed, chrom, start)
    pth <- out_path("intron_orientation", "bed")
    write_tsv(intron_bed, pth, col.names = FALSE)
    built_strandcat$intron_orientation <- pth
    cat(sprintf("intron_orientation: %d intervals (+ %d, - %d)\n", nrow(intron_bed),
                sum(intron_bed$name == "+"), sum(intron_bed$name == "-")))
  } else {
    cat("No strand-unambiguous introns; the intron_orientation annotation is skipped.\n")
  }

  if (has_file("genes6.bed")) {
    genes_named <- fread(work_file("genes6.bed"), header = FALSE, showProgress = FALSE,
                         col.names = c("chrom", "start", "end", "gene", "gene_id", "gene_strand"))
    genes_named[, chrom := as.character(chrom)]
    write_tsv(genes_named, out_path("genes_named", "bed"), col.names = FALSE)
  }

  if (has_file("links.txt")) {
    links <- fread(work_file("links.txt"), header = FALSE, showProgress = FALSE,
                   col.names = c("chrom", "start", "end", "sva_id",
                                 "g_chrom", "g_start", "g_end", "gene", "gene_id", "gene_strand")) %>%
      as_tibble() %>%
      mutate(distance = pmax(0L, pmax(g_start - end, start - g_end))) %>%
      distinct(sva_id, gene, gene_id, gene_strand, distance) %>%
      as.data.table()
    write_tsv(links, out_path("sva_gene_links"))
    cat(sprintf("sva_gene_links: %d links within %s bp (%d SVAs, %d genes)\n", nrow(links), params$link_window,
                n_distinct(links$sva_id), n_distinct(links$gene)))
  }
}

## ==== 7. Replication timing and fork directionality =================================
## Replication timing. There is no hg38 Repli-seq track to download (UCSC has
## none; every ENCODE replication-timing file is hg19 or mm9), so unless
## --repliseq-bigwig is given it is built from 4D Nucleome 2-stage Repli-seq
## (natively GRCh38) as RT = log2(early / late) over the 4DN 5 kb bins.
build_repliseq <- function(early_path, late_path, eps = 1e-3) {
  early <- read_bedgraph(early_path) %>% rename(early = value)
  late <- read_bedgraph(late_path) %>% rename(late = value)
  inner_join(early, late, by = c("chrom", "start", "end")) %>%
    filter(early + late > 0) %>%
    mutate(value = log2((early + eps) / (late + eps))) %>%
    arrange(chrom, start) %>%
    select(chrom, start, end, value)
}

repliseq_path <- fetch(params$repliseq_bigwig, "repliseq")
if (is.null(repliseq_path)) {
  ep <- fetch(params$repliseq_early, "repliseq early fraction")
  lp <- fetch(params$repliseq_late, "repliseq late fraction")
  if (!is.null(ep) && !is.null(lp)) {
    rt <- build_repliseq(ep, lp)
    repliseq_path <- out_path("repliseq_RT_hg38", "bedgraph")
    write_tsv(rt, repliseq_path, col.names = FALSE)
    cat(sprintf("Built replication timing from 4DN Repli-seq: %d bins, RT range %.2f to %.2f -> %s\n",
                nrow(rt), min(rt$value), max(rt$value), repliseq_path))
  }
}
if (!is.null(repliseq_path)) IN$numeric$repliseq <- repliseq_path
IN$numeric <- Filter(Negate(is.null), lapply(setNames(names(IN$numeric), names(IN$numeric)),
                                             function(n) as_bedgraph(IN$numeric[[n]], n)))

## Replication fork directionality (OK-seq RFD). No hg38 RFD track exists
## anywhere (Petryk 2016 deposited raw FASTQ only, SRA SRP065949; every
## processed human RFD is hg19). Routes, in order:
##   1. --rfd-bedgraph: an hg38 RFD track you already have;
##   2. --okseq-watson + --okseq-crick: your own strand-split OK-seq coverage,
##      RFD = (C - W) / (C + W) per bin (Petryk 2016), the cleanest route;
##        samtools view -b -F 0x10 okseq.bam | bedtools genomecov -ibam - -bg > watson.bedgraph
##        samtools view -b -f 0x10 okseq.bam | bedtools genomecov -ibam - -bg > crick.bedgraph
##   3. the published HeLa S3 hg19 track (GSM5814025), lifted to hg38 -- the
##      default, but HeLa and lifted, so a coarse covariate only.
## Positive RFD means rightward-moving forks predominate.
build_rfd <- function(watson_path, crick_path, bin_size, min_cov = 10) {
  bin_size <- as.numeric(bin_size)   # integer bins overflow past 2^31
  ## each interval's coverage spread over fixed bins, weighted by overlap
  binned <- function(path, nm) {
    d <- read_bedgraph(path) %>% filter(chrom %in% genome$chrom)
    b0 <- d$start %/% bin_size
    n <- as.integer((d$end - 1) %/% bin_size - b0 + 1)
    idx <- rep.int(seq_len(nrow(d)), n)
    bin <- sequence(n, from = b0)
    tibble(chrom = d$chrom[idx], bin = bin,
           bp = pmin(d$end[idx], (bin + 1) * bin_size) - pmax(d$start[idx], bin * bin_size),
           cov = d$value[idx]) %>%
      group_by(chrom, bin) %>%
      summarise(!!nm := sum(as.numeric(cov) * bp), .groups = "drop")
  }
  full_join(binned(watson_path, "W"), binned(crick_path, "C"), by = c("chrom", "bin")) %>%
    mutate(W = coalesce(W, 0), C = coalesce(C, 0)) %>%
    filter(W + C >= min_cov) %>%
    mutate(value = (C - W) / (C + W), start = bin * bin_size, end = (bin + 1) * bin_size) %>%
    arrange(chrom, start) %>%
    select(chrom, start, end, value)
}

rfd_path <- IN$rfd
rfd_provenance <- if (!is.null(rfd_path)) "supplied hg38 track" else NA_character_
if (is.null(rfd_path)) {
  wp <- fetch(params$okseq_watson, "okseq_watson")
  cp <- fetch(params$okseq_crick, "okseq_crick")
  if (!is.null(wp) && !is.null(cp)) {
    rfd_tab <- build_rfd(wp, cp, p_int(params$rfd_bin))
    rfd_path <- out_path("okseq_RFD_hg38", "bedgraph")
    write_tsv(rfd_tab, rfd_path, col.names = FALSE)
    rfd_provenance <- "built from your OK-seq coverage (native hg38)"
    cat(sprintf("Built RFD from OK-seq coverage: %d bins of %s bp, median |RFD| = %.3f -> %s\n",
                nrow(rfd_tab), params$rfd_bin, median(abs(rfd_tab$value)), rfd_path))
  } else {
    lifted <- liftover_bedgraph(as_bedgraph(fetch(params$rfd_hg19_bigwig, "published hg19 RFD track"), "hg19 RFD"),
                                fetch(params$liftover_chain, "hg19->hg38 chain"), "hg19 RFD")
    if (!is.null(lifted)) {
      rfd_tab <- read_bedgraph(lifted) %>% filter(chrom %in% genome$chrom)
      rfd_path <- out_path("okseq_RFD_hg38", "bedgraph")
      write_tsv(rfd_tab, rfd_path, col.names = FALSE)
      rfd_provenance <- "published HeLa S3 hg19 RFD, lifted to hg38"
      cat(sprintf(paste("RFD: lifted the published HeLa S3 hg19 track, %d intervals -> %s.",
                        "HeLa and lifted, so use it as a coarse covariate; prefer --okseq-watson/--okseq-crick.\n"),
                  nrow(rfd_tab), rfd_path))
    } else {
      cat("No RFD track could be obtained, so there is no leading/lagging annotation.",
          "Give --rfd-bedgraph, or --okseq-watson and --okseq-crick.\n")
    }
  }
}
if (!is.na(rfd_provenance)) cat("RFD provenance:", rfd_provenance, "\n")

## ==== 8. Annotations ================================================================
LOG("building annotations")
anns <- list()

## One failing source costs that annotation only.
try_build <- function(what, expr) {
  tryCatch(expr, error = function(e) {
    cat(sprintf("WARNING: could not build %s: %s -- skipped.\n", what, conditionMessage(e)))
    NULL
  })
}

## GC per bin (numerator G+C, denominator A+C+G+T), from a GC bedGraph or from
## the reference via bedtools nuc. UCSC gc5Base stores GC as a PERCENT, so a
## track with values above 1.5 is rescaled to a fraction.
gc_bins_from_bedgraph <- function(path) {
  bg <- read_bedgraph(path)[chrom %in% genome$chrom]
  if (max(bg$value, na.rm = TRUE) > 1.5) {
    bg[, value := value / 100]
    cat("GC track looked like percent (max > 1.5), rescaled to a fraction.\n")
  }
  bg <- to_global(bg, genome)
  data.table(gs = bg$gs, num = bg$value * (bg$ge - bg$gs), den = as.numeric(bg$ge - bg$gs))
}
gc_bins_from_reference <- function(fasta) {
  bt <- which_tool("bedtools")
  if (is.null(bt)) stop("GC needs --gc-bedgraph, or --reference with bedtools on PATH")
  ## bedtools cannot index a gzip FASTA (UCSC ships gzip, not bgzip): decompress
  ## once into the cache. hg38 is ~3 GB uncompressed.
  if (grepl("\\.gz$", fasta)) {
    plain <- file.path(CACHE, sub("\\.gz$", "", basename(fasta)))
    if (!file.exists(plain) || file.size(plain) == 0) {
      LOG("decompressing the reference FASTA")
      run_shell(sprintf("gzip -dc %s", shQuote(fasta)), plain)
    }
    fasta <- plain
  }
  gfile <- work_file("gc_genome.txt"); wfile <- work_file("gc_windows.bed"); nfile <- work_file("gc_nuc.txt")
  write_tsv(genome[, .(chrom, size = sprintf("%.0f", size))], gfile, col.names = FALSE)
  run_tool(bt, c("makewindows", "-g", gfile, "-w", p_int(params$gc_bin_size)), stdout = wfile)
  run_tool(bt, c("nuc", "-fi", fasta, "-bed", wfile), stdout = nfile)
  nu <- fread(nfile, showProgress = FALSE)
  setnames(nu, 1:3, c("chrom", "start", "end"))
  nu[, chrom := as.character(chrom)]
  count_col <- function(base) {
    cn <- names(nu)
    nu[[cn[grepl(paste0("num_", base, "_count$"), cn) | grepl(paste0("^\\d+_num_", base), cn)][1]]]
  }
  A <- count_col("A"); C <- count_col("C"); G <- count_col("G"); T_ <- count_col("T")
  nu <- to_global(nu, genome)
  data.table(gs = nu$gs, num = as.numeric(C + G), den = as.numeric(A + C + G + T_))
}
gc_bins <- (if (!is.null(IN$gc_bedgraph)) gc_bins_from_bedgraph(IN$gc_bedgraph)
             else if (!is.null(IN$reference)) gc_bins_from_reference(IN$reference)
             else NULL)
if (!is.null(gc_bins)) {
  anns$gc <- try_build("the GC annotation",
                       make_gc_annotation("gc", gc_bins, genome, window = p_int(params$gc_window),
                                          bin_size = p_int(params$gc_bin_size)))
}

## Gene context, constraint and deserts. Levels resolve BY PRIORITY (first
## match wins), equivalent to the notebook's summarize(annot = min(annot)).
for (nm in names(built_cat)) {
  spec <- built_cat[[nm]]
  bed <- read_bed(spec$path)
  levels_present <- intersect(spec$levels, unique(bed$name))
  if (!length(levels_present)) { cat(sprintf("WARNING: %s has none of its levels; skipped.\n", nm)); next }
  anns[[nm]] <- make_category_annotation(nm, bed, genome, levels_present, spec$background)
}
if (!is.null(anns$gene_context) && !is.null(anns$gene_constraint)) {
  ck <- check_constraint_split(anns$gene_context, anns$gene_constraint, S$g, chrom_index(genome, S$g), S$half)
  cat(sprintf("gene_constraint vs gene_context on the %d SVAs: %d disagree%s\n", ck$n, ck$mismatches,
              if (ck$mismatches) " -- the constraint split is NOT a subdivision of gene_context; do not interpret it"
              else " (an exact subdivision, as required)"))
}

## Intron orientation: the level is the host gene's strand, compared at lookup
## time with the element's own. Gated on gene_context's resolved `intron`, so
## sense + antisense add up to gene_context = intron exactly.
for (nm in names(built_strandcat)) {
  bed <- read_bed(built_strandcat[[nm]])
  gate <- if (!is.null(anns$gene_context)) list(ann = anns$gene_context, level = "intron") else NULL
  if (is.null(gate)) {
    cat(sprintf("NOTE: %s has no gene_context to gate on, so wide elements can exceed the intronic count.\n", nm))
  }
  anns[[nm]] <- try_build(paste("the", nm, "annotation"), make_strand_category_annotation(nm, bed, genome, gate = gate))
}

## Extra categorical BEDs. --category-levels gives each one's priority order,
## pipe-separated ("my_context=CDS|UTR|exon|intron"); without it levels are
## alphabetical, rarely the right priority for nested features.
cat_levels <- p_map(params$category_levels)
cat_backgrounds <- p_map(params$category_backgrounds)
for (nm in names(IN$category)) {
  bed <- read_bed(IN$category[[nm]])
  lv <- if (!is.null(cat_levels[[nm]])) p_vec(cat_levels[[nm]], "|") else sort(unique(bed$name))
  levels_present <- intersect(lv, unique(bed$name))
  if (!length(levels_present)) { cat(sprintf("WARNING: %s has none of the requested levels; skipped.\n", nm)); next }
  bg <- if (!is.null(cat_backgrounds[[nm]])) cat_backgrounds[[nm]] else paste0("no_", nm)
  anns[[nm]] <- make_category_annotation(nm, bed, genome, levels_present, bg)
}

## Numeric tracks: replication timing and any --numeric-annotations.
for (nm in names(IN$numeric)) {
  anns[[nm]] <- make_numeric_annotation(nm, read_bedgraph(IN$numeric[[nm]])[chrom %in% genome$chrom], genome)
}

## Leading/lagging strand: the RFD track plus each element's orientation.
if (!is.null(rfd_path)) {
  if (all(is.na(S$strand))) {
    cat("WARNING: an RFD track is available but no SVA has an orientation; leading/lagging skipped.\n")
  } else {
    anns$lagging_strand <- make_strand_match_annotation(
      "lagging_strand", read_bedgraph(rfd_path)[chrom %in% genome$chrom], genome,
      flip = p_bool(params$rfd_flip), min_abs = p_num(params$rfd_min_abs))
    cat(sprintf("Leading/lagging: |RFD| < %s is ambiguous, flip = %s.\n", params$rfd_min_abs, params$rfd_flip))
  }
}

if (!length(anns)) {
  stop("No annotations were built. Give at least one of --reference (or --gc-bedgraph), --gencode-gtf, ",
       "--numeric-annotations, --rfd-bedgraph or --category-annotations.")
}
cat(sprintf("Annotations: %s\n", paste(sprintf("%s (%s)", names(anns), vapply(anns, `[[`, "", "type")), collapse = ", ")))

## Is the RFD track the right way round? Sign conventions differ between
## OK-seq pipelines, and an inverted track silently swaps lagging and leading.
## RFD rises across initiation zones, which sit in early-replicating DNA, so
## the local slope of RFD should correlate POSITIVELY with replication timing.
## Deliberately not checked on the SVAs themselves: that is the biology under test.
rfd_polarity <- NULL
if (!is.null(rfd_path) && !is.null(anns$repliseq) && !is.null(anns$lagging_strand)) {
  slopes <- read_bedgraph(rfd_path) %>%
    filter(chrom %in% genome$chrom) %>%
    arrange(chrom, start) %>%
    group_by(chrom) %>%
    mutate(prev_end = lag(end), prev_val = lag(value)) %>%
    ungroup() %>%
    filter(!is.na(prev_val), start - prev_end < 2 * (end - start)) %>%   # adjacent bins only
    mutate(slope = value - prev_val)
  if (nrow(slopes) > 100) {
    mid <- gpos(genome, slopes$chrom, (slopes$start + slopes$end) %/% 2L)
    d <- tibble(slope = slopes$slope, rt = lookup_annotation(anns$repliseq, mid, chrom_index(genome, mid))) %>%
      filter(is.finite(slope), is.finite(rt))
    if (nrow(d) > 100) {
      ct <- suppressWarnings(cor.test(d$slope, d$rt, method = "spearman"))
      rho <- unname(ct$estimate)
      verdict <- case_when(is.na(rho) ~ "INCONCLUSIVE",
                           rho > 0.02 ~ "consistent with the assumed polarity",
                           rho < -0.02 ~ "LOOKS INVERTED -- set --rfd-flip true",
                           TRUE ~ "INCONCLUSIVE (correlation near zero)")
      rfd_polarity <- data.table(n_bins = nrow(d), spearman_slope_vs_timing = rho, p = ct$p.value,
                                 rfd_flip_setting = params$rfd_flip, verdict = verdict)
      show_tab(rfd_polarity, caption = "RFD polarity check (slope of RFD vs replication timing should be positive)")
      write_tsv(rfd_polarity, out_path("rfd_polarity_check"))
    }
  }
  if (is.null(rfd_polarity)) cat("RFD polarity could not be checked (too few adjacent bins); verify the sign by hand.\n")
}

## ==== 9. Allele counts and frequencies ==============================================
## As in variant_overview.Rmd: AC = hap1 + hap2 summed over individuals per
## clustered SVA, chrX genotypes missing for individuals whose sex is not XX or
## XY; AN = 2 per individual on autosomes, 2 per XX and 1 per XY on chrX, 1 per
## XY on chrY. Without metadata only the autosomal AN is known.
LOG("allele counts")
is_carried <- function(x) as.numeric(toupper(stri_trim_both(x)) %in% c("TRUE", "T", "1"))
read_genotypes <- function(path, id_fmt) {
  as_tibble(fread_any(path, sep = "\t", colClasses = "character", showProgress = FALSE)) %>%
    mutate(sva_id = sprintf(id_fmt, as.integer(cluster_id)))
}
gt <- read_genotypes(IN$genotype_table, "SVA_MEI_cl%05d")
md <- NULL
if (is.null(IN$sample_metadata)) {
  cat("No --sample-metadata: chrX genotypes are not sex-filtered and chrX/chrY AF is missing.\n")
} else {
  ## participant 1404227 has two conflicting sex entries
  md <- as_tibble(fread_any(IN$sample_metadata, sep = "\t", colClasses = "character", showProgress = FALSE)) %>%
    distinct(indiv = participant, sex) %>%
    filter(indiv %in% gt$indiv, indiv != "1404227")
}
sex_exclude <- if (is.null(md)) character() else md$indiv[!md$sex %in% c("XX", "XY")]
n_ind <- if (is.null(md)) n_distinct(gt$indiv) else n_distinct(md$indiv)
AN_AUTO <- 2 * n_ind
AN_X <- if (is.null(md)) NA_real_ else 2 * n_distinct(md$indiv[md$sex == "XX"]) + n_distinct(md$indiv[md$sex == "XY"])
AN_Y <- if (is.null(md)) NA_real_ else as.numeric(n_distinct(md$indiv[md$sex == "XY"]))

## AC / AN / AF per element of one family's genotype table
count_alleles <- function(gt) {
  gt %>%
    mutate(n_carried = if_else(CHROM == "chrX" & indiv %in% sex_exclude, NA_real_, is_carried(hap1) + is_carried(hap2))) %>%
    group_by(CHROM, sva_id) %>%
    summarise(AC = sum(n_carried, na.rm = TRUE), .groups = "drop") %>%
    mutate(AN = case_when(CHROM %in% c("chrX", "X") ~ AN_X, CHROM %in% c("chrY", "Y") ~ AN_Y, TRUE ~ AN_AUTO),
           AF = AC / AN) %>%
    distinct(sva_id, .keep_all = TRUE) %>%
    select(sva_id, AC, AN, AF)
}
allele_counts <- count_alleles(gt)
S <- S %>%
  left_join(allele_counts, by = c("name" = "sva_id")) %>%
  mutate(across(c(AC, AN, AF), ~ if_else(class == PL, as.numeric(.x), NA_real_)))
cat(sprintf("Genotypes: %d SVAs, %d individuals (AN autosomes %d, chrX %s); %d of %d %s SVAs have AC > 0, %d are singletons.\n",
            nrow(allele_counts), n_ind, as.integer(AN_AUTO), format(AN_X), sum(S$AC > 0, na.rm = TRUE),
            sum(S$class == PL), PL, sum(S$AC == 1, na.rm = TRUE)))
if (nrow(OTHER)) {
  ## the same individuals and allele numbers as the SVAs: one cohort, one call set
  other_ac <- bind_rows(lapply(names(OTHER_FAMILIES), function(fam) {
    count_alleles(read_genotypes(OTHER_FAMILIES[[fam]]$genotype_table, OTHER_FAMILIES[[fam]]$id_fmt))
  }))
  OTHER <- OTHER %>% left_join(other_ac, by = c("name" = "sva_id"))
  show_tab(OTHER %>% group_by(family) %>%
             summarise(n = n(), with_ac = sum(AC > 0, na.rm = TRUE), singletons = sum(AC == 1, na.rm = TRUE)),
           caption = "Alu / L1 allele counts")
}

## ==== 10. Structural class from the MEI structural catalog ===========================
## The catalog's variant_id is the MEI table's var_id by default, or
## "<MEI_type>_cl<cluster_id>" when it was run with MeiVariantCol = cluster_id.
## Each catalog variant is mapped to its clustered SVA by cluster id, then by the
## MEI bed's ID, then by position; an SVA spanning several catalog variants takes
## the modal structure and majority flags, each variant weighted by the
## haplotypes carrying it (as the catalog rolls sequences up). structure_group
## takes transduction first, because a transduced element can also be
## full-length or truncated. catalog_subfamily is the weighted modal catalog
## subfamily: the catalog recognises SVA_F1 by its MAST2 exon 1, which the MEI
## bed's labels (RepeatMasker has no F1 consensus) call SVA_F or SVA_E.
catalog_structure <- function(path, mei_bed, catalog_families = "SVA", id_fmt = "SVA_MEI_cl%05d") {
  flag_cols <- c("full_length", "five_prime_truncated", "three_prime_truncated",
                 "three_prime_transduction", "five_prime_transduction")
  pv <- as_tibble(fread_any(path, sep = "\t", showProgress = FALSE))
  miss <- setdiff(c("variant_id", "structure", flag_cols), names(pv))
  if (length(miss)) stop(path, ": not a structure_per_variant table, missing ", paste(miss, collapse = ", "))
  if ("family" %in% names(pv)) pv <- pv %>% filter(toupper(family) %in% catalog_families)
  if (!"allele_count" %in% names(pv)) pv$allele_count <- 1
  if (!"subfamily" %in% names(pv)) pv$subfamily <- NA_character_
  pv <- pv %>% mutate(variant_id = as.character(variant_id), across(all_of(flag_cols), as.logical))

  b <- fread_any(mei_bed, sep = "\t", header = FALSE, showProgress = FALSE)
  setnames(b, seq_along(MEI_COLS), MEI_COLS)
  bed <- as_tibble(b) %>%
    transmute(ID = as.character(ID), chrom = as.character(chrom), start = as.integer(start),
              cluster_id = as.integer(cluster_id), MEI_type = as.character(MEI_type)) %>%
    filter(!is.na(cluster_id))
  ## with MeiVariantCol = cluster_id the catalog names a variant <MEI_type>_cl<id>
  name_re <- paste0("(?i)^(", paste(unique(c(bed$MEI_type, catalog_families)), collapse = "|"), ")_cl0*(\\d+)$")
  by_id <- bed %>% distinct(ID, .keep_all = TRUE) %>% select(variant_id = ID, cluster_by_id = cluster_id)
  ## a catalog `pos` is the bed start (0-based) or start + 1; only unambiguous sites
  by_pos <- bind_rows(bed %>% transmute(chrom, pos = start, cluster_id),
                      bed %>% transmute(chrom, pos = start + 1L, cluster_id)) %>%
    distinct() %>%
    group_by(chrom, pos) %>%
    filter(n_distinct(cluster_id) == 1) %>%
    slice(1) %>%
    ungroup() %>%
    rename(cluster_by_pos = cluster_id)

  pv <- pv %>%
    mutate(cluster_by_name = as.integer(stri_match_first_regex(variant_id, name_re)[, 3])) %>%
    left_join(by_id, by = "variant_id")
  pv <- if (all(c("chrom", "pos") %in% names(pv))) {
    pv %>% mutate(chrom = as.character(chrom), pos = as.integer(pos)) %>% left_join(by_pos, by = c("chrom", "pos"))
  } else {
    pv %>% mutate(cluster_by_pos = NA_integer_)
  }
  pv <- pv %>%
    mutate(cluster_id = coalesce(cluster_by_name, cluster_by_id, cluster_by_pos),
           matched_by = case_when(!is.na(cluster_by_name) ~ "variant_id is a cluster id",
                                  !is.na(cluster_by_id) ~ "variant_id = MEI bed ID",
                                  !is.na(cluster_by_pos) ~ "chrom + position",
                                  TRUE ~ "unmatched"))
  show_tab(count(pv, matched_by), caption = sprintf("Catalog %s variants matched to clustered insertions",
                                                     paste(catalog_families, collapse = "/")))

  majority <- function(x, w) sum(w * (x %in% TRUE)) / sum(w) >= 0.5
  pv %>%
    filter(!is.na(cluster_id)) %>%
    mutate(sva_id = sprintf(id_fmt, cluster_id), w = coalesce(as.numeric(allele_count), 0)) %>%
    group_by(sva_id) %>%
    mutate(w = if (sum(w) > 0) w else rep(1, n())) %>%   # nothing carried: equal weights
    summarise(structure = names(which.max(tapply(w, coalesce(structure, "NA"), sum))),
              catalog_subfamily = names(which.max(tapply(w, coalesce(subfamily, "NA"), sum))),
              is_full_length = majority(full_length, w),
              is_truncated = majority(five_prime_truncated | three_prime_truncated, w),
              has_transduction = majority(three_prime_transduction | five_prime_transduction, w),
              n_catalog_variants = n(),
              .groups = "drop") %>%
    mutate(structure_group = case_when(has_transduction ~ "transduction",
                                       is_full_length ~ "full_length",
                                       is_truncated ~ "truncated",
                                       TRUE ~ "other"))
}

struct <- NULL
if (!is.null(IN$catalog)) {
  LOG("structural classes from the catalog")
  struct <- catalog_structure(IN$catalog, IN$mei_bed)
  cat(sprintf("Structural class for %d of %d %s SVAs:\n",
              sum(S$class == PL & S$name %in% struct$sva_id), sum(S$class == PL), PL))
  show_tab(struct %>% filter(sva_id %in% S$name) %>% count(structure_group, sort = TRUE))
  cat(sprintf("%d %s SVAs relabelled SVA_F1 from the catalog (MAST2 exon 1 at the 5' end).\n",
              sum(S$name %in% struct$sva_id[struct$catalog_subfamily %in% "SVA_F1"]), PL))
  if (nrow(OTHER)) {
    other_struct <- bind_rows(lapply(names(OTHER_FAMILIES), function(fam) {
      f <- OTHER_FAMILIES[[fam]]
      catalog_structure(IN$catalog, f$bed, f$catalog_families, f$id_fmt)
    }))
    show_tab(OTHER %>% left_join(other_struct, by = c("name" = "sva_id")) %>% count(family, structure_group),
             caption = "Alu / L1 structural classes")
  }
} else {
  cat("No --catalog-per-variant: the structure columns are left out.\n")
}

## ==== 11. Outputs ===================================================================
LOG("per-SVA table")
## Every annotation's value per SVA, labelled as the notebook counts it. A name
## that clashes with another column gets an `ann_` prefix; the mapping travels
## in the bundle.
ann_values <- as_tibble(element_annotations(anns, S$g, genome, S$strand, S$half))
ann_names <- names(ann_values)
STRUCT_COLS <- c("structure", "is_full_length", "is_truncated", "has_transduction", "n_catalog_variants", "structure_group",
                 "catalog_subfamily")
LENGTH_COLS <- intersect(c("n_obs", "n_indiv", "n_alleles", "ins_length", "n_hexamer", "n_VNTR", LENGTH_STATS), names(S))
LINK_COLS <- c("n_linked_genes", "linked_genes", "nearest_linked_gene_distance")
clash <- intersect(ann_names, c("sva_id", "chrom", "start", "end", "strand", "subfamily", "subfamily_mei", "class", "is_polymorphic",
                                "is_singleton", "ac", "an", "af", "width", STRUCT_COLS, LENGTH_COLS, LINK_COLS))
names(ann_values) <- ifelse(ann_names %in% clash, paste0("ann_", ann_names), ann_names)

tab <- S %>%
  transmute(sva_id = name, chrom, start, end, strand, subfamily, class,
            is_polymorphic = class == PL,
            is_singleton = class == PL & AC %in% 1,
            ac = AC, an = AN, af = AF, width,
            across(all_of(LENGTH_COLS))) %>%
  bind_cols(ann_values)
if (!is.null(struct)) {
  tab <- tab %>%
    left_join(struct, by = "sva_id") %>%
    mutate(subfamily_mei = subfamily,
           subfamily = if_else(catalog_subfamily %in% "SVA_F1", "SVA_F1", subfamily)) %>%
    relocate(subfamily_mei, .after = subfamily) %>%
    relocate(all_of(STRUCT_COLS), .after = width)
}
if (!is.null(links) && nrow(links)) {
  linked <- as_tibble(links) %>%
    arrange(distance, gene) %>%
    group_by(sva_id) %>%
    summarise(n_linked_genes = n_distinct(gene),
              linked_genes = paste(unique(gene), collapse = ";"),
              nearest_linked_gene_distance = min(distance),
              .groups = "drop")
  tab <- tab %>%
    left_join(linked, by = "sva_id") %>%
    mutate(n_linked_genes = coalesce(n_linked_genes, 0L))
}
tab <- tab %>% arrange(chrom, start)
table_path <- out_path("sva_annotation_table", "tsv.gz")
write_tsv(tab, table_path)

## One table per other family (alu_annotation_table, l1_annotation_table): the
## same per-element columns as the SVA table, minus the SVA repeat lengths.
if (nrow(OTHER)) {
  other_ann <- as_tibble(element_annotations(anns, OTHER$g, genome, OTHER$strand, OTHER$half))
  names(other_ann) <- names(ann_values)
  other_tab <- OTHER %>%
    transmute(mei_id = name, family, chrom, start, end, strand, subfamily, class,
              is_polymorphic = TRUE, is_singleton = AC %in% 1, ac = AC, an = AN, af = AF, width,
              ins_length, n_alleles) %>%
    bind_cols(other_ann)
  if (!is.null(struct)) {
    other_tab <- other_tab %>%
      left_join(other_struct, by = c("mei_id" = "sva_id")) %>%
      relocate(any_of(STRUCT_COLS), .after = width)
  }
  if (!is.null(links) && nrow(links)) {
    other_tab <- other_tab %>%
      left_join(linked, by = c("mei_id" = "sva_id")) %>%
      mutate(n_linked_genes = coalesce(n_linked_genes, 0L))
  }
  for (fam in unique(other_tab$family)) {
    fam_path <- out_path(paste0(tolower(fam), "_annotation_table"), "tsv.gz")
    fam_tab <- other_tab %>% filter(family == fam) %>% arrange(chrom, start)
    write_tsv(fam_tab, fam_path)
    cat(sprintf("Wrote %s (%d %s insertions).\n", fam_path, nrow(fam_tab), fam))
  }
}

## Everything the notebook needs besides the table. Tables are data.tables
## because the notebook indexes them that way.
bundle_path <- out_path("sva_annotations", "rds")
saveRDS(list(anns = anns, genome = genome, excl = excl, incl = incl,
             links = links, genes_named = genes_named,
             funnel = as.data.table(FUNNEL %>% select(stage, class, n)),
             rfd_polarity = rfd_polarity, rfd_provenance = rfd_provenance,
             annotation_columns = setNames(names(ann_values), ann_names),
             lib = LIB, annotation_checks = ANNOTATION_CHECKS,
             params = params, created = Sys.time()),
        bundle_path)
cat(sprintf("Wrote %s (%d SVAs x %d columns) and %s.\n", table_path, nrow(tab), ncol(tab), bundle_path))
show_tab(tibble(column = names(tab), non_missing = vapply(tab, function(v) sum(!is.na(v)), integer(1))),
         caption = "Columns")
LOG("done")
