# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Bioinformatics workflows for genome-wide association work on mobile element insertions (MEIs), specifically SVA elements, in long-read sequencing data. There are two independent WDL workflows, registered for Dockstore, that run on Cromwell/Terra.

## Architecture

Two unrelated pipelines share this repo; they do not call each other.

### Overall Organization
- workflows/ contains any WDLs that will be used for pipelines, subfolders and utility WDL files are allowed
- scripts/ contains any scripts that will be used within custom environments for pipelines, every dockerfile should copy the entire scripts/ folder 
- envs/ contain Dockerfiles for custom environments used in workflows, named flat as `Dockerfile.<env-name>` (not one Dockerfile per subfolder)

All envs should have their own rule in .github/workflows so that they will be updated whenever the relevant dockerfile or any of scripts/* is changed
All workflows should be registered in .dockstore.yml

### 1. ExtractVariants (`workflows/extract_variants.wdl` + `scripts/extract_snps.py`)

Pulls genotypes for a set of variants and samples out of a Hail VDS stored in GCS.

- The WDL has a single task that invokes `/extract_snps.py` inside the Hail Docker image.
- `extract_snps.py` runs on Spark via Hail. Key flow:
  - Reads a **variant list** (one `chrom:pos_ref_alt` string per line) and a **sample list** (one sample ID per line).
  - Bins variants into ≤2 Mb intervals to make `hl.vds.filter_intervals` efficient, then applies `filter_variants` and `filter_samples`.
  - Reconstructs dense genotypes from the VDS local-allele representation: `local_to_global` for AD, `lgt_to_gt` for GT, then `to_dense_mt`.
  - Exports `GT` to `genotypes.tsv` (the workflow's only output).
- The variant-string parsing format (`chrom:pos_ref_alt`) is load-bearing — it is parsed both in Python-native code and again via Hail expressions. Keep both parsers in sync if the format changes.
- Docker image is built from **this repo** (`envs/Dockerfile.hail`, based on `hailgenetics/hail`) and published to Docker Hub as `<DOCKERHUB_USERNAME>/mei-lr-association-hail` (CI names each image `<repo-name>-<env-name>`). The WDL selects the tag via the `ImageTag` input (defaults to `latest`; pass a 7-char commit SHA to pin a specific build).

### 2. SVA_Typer (`workflows/sva_typer.wdl`)

Types the internal repeat structure of SVA elements from a FASTA.

- Runs an external `sva_typer` binary (`--write-query-seq-state`) to emit per-base repeat state, then post-processes with `process_output.py` to get repeat lengths.
- Both the `sva_typer` tool and `process_output.py` live **inside the Docker image** (`ayenkin1871/sva_typer:latest`), NOT in this repo. This repo only contains the WDL wrapper.

### 3. ExtractInsertionMethylation (`workflows/extract_insertion_methylation.wdl` + `scripts/extract_insertion_methylation.py`)

Extracts per-read methylation over insertion sites from a haplotagged PacBio HiFi modBAM.

- Inputs: one indexed modBAM (MM/ML tags + `HP` haplotag) and a BED of insertion loci — the MEI bed layout (see section 4 below): `chrom, start, end, var_id, score, strand, min_len, max_len` (8 columns). `score`/`strand` are present but unused by this workflow; only `min_len`/`max_len` (columns 7/8) bound the expected insertion length in bp. This workflow only makes sense for non-reference (MEI) insertion loci — a reference SVA bed has no `min_len`/`max_len` and isn't a valid `LociBed` here, since there's no insertion CIGAR op to find for something already in the reference. Output: one TSV row per (locus, read).
- `extract_insertion_methylation.py` (pysam) walks each read's CIGAR to find the large `I` (insertion) op anchored at the locus, extracts the inserted bases, and intersects them with the read's `modified_bases` (5mC `m` / 5hmC `h`). Emits the insertion sequence, haplotype, a per-base methylation string, and per-mod summary stats.
- Accepted insertion length is per-locus: an `I` op of length `L` qualifies when `min_len - LenTolerance <= L <= max_len + LenTolerance` (`--len-tolerance`, default 100 bp), replacing the old global `--min-ins-len`/`--max-ins-len`.
- The insertion span itself is a `query_sequence` slice — no reference coordinates apply to inserted bases.
- Key correctness invariant: `read.modified_bases` positions and CIGAR query offsets are both in `query_sequence` coordinates (forward-reference oriented), so they intersect directly — use `modified_bases`, **not** `modified_bases_forward`. Probability = `(qual + 0.5) / 256`; `qual == -1` is a no-call.
- Docker image built from `envs/Dockerfile.python_general` (a general-purpose Python env — python + pysam + pandas + biopython — also used by other pure-Python scripts like `combine_fastas.py`), published to Docker Hub as `<DOCKERHUB_USERNAME>/mei-lr-association-python_general`; WDL selects the tag via the `ImageTag` input (defaults to `latest`). The WDL symlinks `BamIndex` next to `Bam` so pysam finds the `.bai`.

### 4. SVAMEIMethylationWindows (`workflows/sva_mei_methylation_windows.wdl`)

End-to-end pipeline, for one individual, from a reference-SVA bed + an MEI bed to a 100bp-binned methylation-by-distance-from-variant summary, built from several composable sub-workflows:

- The reference SVA bed and MEI bed are two distinct catalogs with different (but overlapping) column layouts. Both carry the variant ID in column 4 and strand in column 6, by convention:
  - SVA bed: `#chrom, start, end, sva_id, score, strand, subfamily` (7 columns).
  - MEI bed: `#chrom, start, end, var_id, score, strand, min_length, max_length` (8 columns).
  - `IndividualVariantBedFromLongGenotypes` (`workflows/individual_variant_bed_from_long_genotypes.wdl`) filters either catalog down to one individual's variants by `var_id` (via a sparse long-form genotype table); it emits matching rows unmodified (awk's bare `print` of the whole line), so it preserves whatever columns the input bed has in their normal positions — including `score`/`strand` — not just the column it matched on.
  - `ExtractInsertionMethylation` (section 3 above) consumes the MEI bed layout specifically (only makes sense for non-reference insertions).
- **CombineAndWindowVariants** (`workflows/combine_and_window_variants.wdl`): concatenates the SVA and MEI beds (selecting columns 1/2/3/4/6 by position from each, so it doesn't care that the two schemas diverge past column 6), then flanks every variant by `FlankSize` (default 10000 bp) on each side via `bedtools slop -g <ChromSizes>` (clamped to chromosome bounds). The pre-flank coordinates and strand are carried through as extra trailing columns (`var_start`, `var_end`) rather than being overwritten, since `bedtools slop` only touches the start/end columns — this is what distance-from-variant is measured against downstream. Output layout: `chrom, start, end, var_id, strand, var_start, var_end` (7 columns).
- **BedtoolsIntersect** (`workflows/bedtools_intersect.wdl`, pre-existing, unmodified) is called once per haplotype — `Hap1MethylBed`/`Hap2MethylBed`, each a haplotype-specific pb-CpG-tools pileup bed (`chrom, start, end, modification_probability, haplotype, coverage, est_mod_site_count, est_unmod_site_count, discrete_modification_score`, 9 columns) — against the windowed variant bed, with `-wa -wb` so both sides' columns are kept.
- **ProcessMethylationWindow** (`workflows/process_methylation_window.wdl` + `scripts/process_methylation_window.R`) is then called once per haplotype, independently, directly on each haplotype's own intersect bed — it is deliberately *not* given a pre-concatenated, pre-tagged two-haplotype bed, so it stays decoupled from how many haplotype files there are (it can equally be run standalone later on a single combined/unphased methylation bed). It bins each CpG site into a fixed `BinSize`-wide tile (default 100 bp) by signed distance from its variant, oriented by strand so downstream is positive, out to `FlankSize` bp each way; a site falling inside the variant's own (pre-flank) span gets the special `intersect` bin instead of a distance tile. Summarizes `num_cpgs`, `avg_coverage`, and `mean_methylation_score` (from `modification_probability`) per `(var_id, haplotype, dist_bin)` — haplotype grouping comes from the methylation bed's own embedded `haplotype` column (not from which file was passed in); the methylation bed's extra `est_mod_site_count`/`est_unmod_site_count`/`discrete_modification_score` columns are parsed (to keep column positions aligned with the variant-bed columns that follow) but not otherwise used. This replaced an earlier version of the script's distance binning, which used a small number of named, uneven categories (`ups_1-5kb`, etc.) instead of fixed-width tiles.
- The two per-haplotype summary TSVs (each already carrying a header and its own `haplotype` column) are concatenated afterward via `ConcatFiles.ConcatenateDelim` (`workflows/utils/concat_files.wdl`, `HasHeader = true`) into the final merged summary. The per-haplotype summaries are also exposed directly as workflow outputs (`Hap1MethylationWindowSummary`/`Hap2MethylationWindowSummary`), not just the merged one.
- Docker images: `CombineAndWindowVariants`/`BedtoolsIntersect` run in the public `quay.io/biocontainers/bedtools` image (same as the sibling `AoU_mQTL_analysis` repo's binning tasks); `ConcatenateDelim` and `ProcessMethylationWindow` run in this repo's `python_general`/`r_analysis` images respectively.

### nbsync (`tools/nbsync`)

Syncs `notebooks/` in this repo with the `notebooks/` folder in the Terra workspace bucket (the default location Terra puts notebooks in). That bucket folder also holds many other scripts/files unrelated to this repo — nbsync only ever touches files that have been explicitly **added**, so it won't pull in everything else sitting there.

- **Git-like tracking, no commit step**: a file sitting in `notebooks/` is just "untracked" until you run `tools/nbsync add <file> [<file> ...]` (errors if the file isn't present locally yet — use `pull --add` instead for a file that only exists in the bucket so far). Once added, it's immediately eligible for `pull`/`push` — there's no separate commit/staging phase.
- `tools/nbsync status` — like `git status`: lists `Added:` files (each with a one-line sync status: in sync / local changes / remote changes / conflict) and `Untracked:` files (present in `notebooks/` but never added). Read-only except for a one-time silent migration of `.nbsync_state.json` to the current format (see below).
- `tools/nbsync pull` downloads (overwrites local with) the bucket copy of every **added** file.
- `tools/nbsync pull --add <file> [<file> ...]` adds and immediately pulls specific named files from the bucket, even if not already added — this is how you start tracking a new notebook/script that already exists on Terra.
- `tools/nbsync push` uploads every **added** file present in `notebooks/` to the bucket, overwriting the bucket copy.
- Bucket name is read from `.nbsync.secret` at the repo root — a single line with just the bucket name (no `gs://` prefix), gitignored and excluded from Claude's `Read` tool via `.claude/settings.local.json`. Create it yourself (not via Claude) with e.g. `echo "fc-secure-xxxxxxxx" > .nbsync.secret`.
- Uses `gsutil cp`/`gsutil stat` under the hood (one invocation per file — each pays gsutil's own CLI startup cost, so a `pull`/`push`/`status` over many files can take noticeably long); requires `gcloud`/`gsutil` to be authenticated already.
- **State file / migration**: tracking + sync state both live in `.nbsync_state.json` at the repo root (gitignored, per-checkout — not shared across collaborators' clones): `{"added": [...], "files": {name: {"hash", "pulled_at"}}}`. An older flat-format state file (from before the explicit `add` step existed, `{name: {"hash", "pulled_at"}}` with no top-level `"added"`/`"files"`) is migrated automatically and transparently the first time *any* command loads it — every file physically present in `notebooks/` at that moment is treated as already added (so upgrading doesn't silently stop syncing anything), and the migrated structure is written back to disk immediately (not just held in memory) so the union-with-current-`notebooks/`-contents logic only ever runs once; otherwise a stale on-disk old-format file would keep re-adopting newly-dropped, genuinely-untracked files as "added" on every subsequent read.
- **Conflict detection**: every successful pull records the file's content (md5) hash in the state file's `"files"` map. `push` only goes through if the bucket's current hash still matches that recorded hash (nobody else pushed since your last pull); otherwise it's blocked with a "pull before push" message, per file. `pull` only succeeds outright if your local file's hash still matches the recorded hash (you haven't edited it locally since your last pull); if local is edited but remote hasn't moved, the pull is skipped (push to save your edits instead). If *both* local and remote have diverged from the recorded hash, you get an interactive prompt: force-pull (overwrite local), cancel, or vimdiff-merge (opens `vimdiff <local> <downloaded-remote>`; once you save and quit, the recorded hash is pinned to the remote's, so a subsequent `push` goes through even though local content differs from both the old and new remote). A newly-added file with no recorded hash yet is only ever treated as a conflict if local and remote actually differ — if they already agree (e.g. already synced before this hash-tracking existed) it's silently adopted into tracking.
- **Error vs. "not found" distinction**: `gsutil stat` exiting nonzero only means the object is genuinely missing when its stderr is a clean one-line message; a Python traceback (auth expired, network down, etc.) raises `RemoteCheckError` instead of being treated as "file not found in bucket" or silently folded into a hash mismatch — the command aborts immediately with a clear message (e.g. suggesting `gcloud auth login`) rather than misreporting every file as changed/missing and burning minutes retrying the same failure file-by-file. This was a real bug caught in practice: an expired `gcloud` reauth session made every tracked file falsely show "remote has changed" until fixed.

## Common Commands

There is no build/test/lint tooling in this repo — it is WDL + a single Python script + a Dockerfile.

Build the Hail Docker image locally (note: Docker context is repo root, Dockerfile path is explicit):
```
docker build -f envs/Dockerfile.hail -t mei-lr-association .
```

Validate / run WDLs locally (requires womtool/cromwell, not vendored):
```
womtool validate workflows/extract_variants.wdl
```

## CI

Each env has its own workflow (`.github/workflows/docker-image.yml` for hail, `python-general-docker-image.yml` for python_general, `r-analysis-docker-image.yml` for r_analysis, `bioinformatics-docker-image.yml` for bioinformatics) that builds and pushes to Docker Hub on push/PR to `main`/`develop`, **only when `scripts/**` or that env's `Dockerfile.<env-name>` change**. Editing the WDLs alone will not trigger a rebuild. Images are named `<DOCKERHUB_USERNAME>/<repo-name>-<env-name>` and tagged `latest` + the 7-char commit SHA. Auth uses the `DOCKERHUB_USERNAME` repo **variable** and the `DOCKERHUB_TOKEN` **secret**. Because all workflows watch `scripts/**`, a script change rebuilds every image — the `<env-name>` suffix keeps their names distinct.

## Gotchas

- **Dockstore path casing**: `.dockstore.yml` points ExtractVariants at `/workflows/ExtractVariants.wdl`, but the file on disk is `workflows/extract_variants.wdl`. If Dockstore registration fails, reconcile the casing.
- Spark/Hail resource settings are hardcoded in `extract_snps.py` (`hl.init` spark_conf) and must be kept consistent with the WDL `runtime` block (`memory`, `cpu`, `disks`). `SPARK_LOCAL_DIRS`/`spark.local.dir` both point at `/cromwell_root` for shuffle spill.
