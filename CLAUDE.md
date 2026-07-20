# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Bioinformatics workflows for genome-wide association work on mobile element insertions (MEIs), specifically SVA elements, in long-read sequencing data. There are two independent WDL workflows, registered for Dockstore, that run on Cromwell/Terra.

## Architecture

Two unrelated pipelines share this repo; they do not call each other.

### Overall Organization
- workflows/ contains any WDLs that will be used for pipelines, subfolders and utility WDL files are allowed
- scripts/ contains any scripts that will be used within custom environments for pipelines, every dockerfile should copy the entire scripts/ folder 
- envs/ contain Dockerfiles for custom environments used in workflows

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
- Docker image is built from **this repo** (`envs/hail/Dockerfile`, based on `hailgenetics/hail`) and published to `ghcr.io/alyenkin/mei-lr-association`. The WDL selects the tag via the `Branch` input (defaults to `main`).

### 2. SVA_Typer (`workflows/sva_typer.wdl`)

Types the internal repeat structure of SVA elements from a FASTA.

- Runs an external `sva_typer` binary (`--write-query-seq-state`) to emit per-base repeat state, then post-processes with `process_output.py` to get repeat lengths.
- Both the `sva_typer` tool and `process_output.py` live **inside the Docker image** (`ayenkin1871/sva_typer:latest`), NOT in this repo. This repo only contains the WDL wrapper.

### 3. ExtractInsertionMethylation (`workflows/extract_insertion_methylation.wdl` + `scripts/extract_insertion_methylation.py`)

Extracts per-read methylation over insertion sites from a haplotagged PacBio HiFi modBAM.

- Inputs: one indexed modBAM (MM/ML tags + `HP` haplotag) and a BED of insertion loci (~2–4 kb each). Output: one TSV row per (locus, read).
- `extract_insertion_methylation.py` (pysam) walks each read's CIGAR to find the large `I` (insertion) op anchored at the locus, extracts the inserted bases, and intersects them with the read's `modified_bases` (5mC `m` / 5hmC `h`). Emits the insertion sequence, haplotype, a per-base methylation string, and per-mod summary stats.
- Key correctness invariant: `read.modified_bases` positions and CIGAR query offsets are both in `query_sequence` coordinates (forward-reference oriented), so they intersect directly — use `modified_bases`, **not** `modified_bases_forward`. Probability = `(qual + 0.5) / 256`; `qual == -1` is a no-call.
- Docker image built from `envs/methylation/Dockerfile` (python + pysam), published to `ghcr.io/alyenkin/insertion-methylation`; WDL selects the tag via the `Branch` input. The WDL symlinks `BamIndex` next to `Bam` so pysam finds the `.bai`.

## Common Commands

There is no build/test/lint tooling in this repo — it is WDL + a single Python script + a Dockerfile.

Build the Hail Docker image locally (note: Docker context is repo root, Dockerfile path is explicit):
```
docker build -f envs/hail/Dockerfile -t mei-lr-association .
```

Validate / run WDLs locally (requires womtool/cromwell, not vendored):
```
womtool validate workflows/extract_variants.wdl
```

## CI

`.github/workflows/docker-image.yml` builds and pushes the Hail image to `ghcr.io/alyenkin/mei-lr-association` on push/PR to `main` (and pushes to `develop`), **only when files under `scripts/` change**. Editing the Dockerfile or WDLs alone will not trigger a rebuild — touch something in `scripts/` (or trigger manually) to publish a new image.

## Gotchas

- **Dockstore path casing**: `.dockstore.yml` points ExtractVariants at `/workflows/ExtractVariants.wdl`, but the file on disk is `workflows/extract_variants.wdl`. If Dockstore registration fails, reconcile the casing.
- Spark/Hail resource settings are hardcoded in `extract_snps.py` (`hl.init` spark_conf) and must be kept consistent with the WDL `runtime` block (`memory`, `cpu`, `disks`). `SPARK_LOCAL_DIRS`/`spark.local.dir` both point at `/cromwell_root` for shuffle spill.
