version 1.0

## Computes per-base coverage across an assembly-to-reference alignment BAM
## (e.g. contigs aligned to the reference with minimap2 asm5), reporting two
## depths side by side so gaps vs. deletions can be told apart downstream:
##   - `depth`:           count of reads with an actual aligned base at this
##                        position (excludes reads whose CIGAR has a deletion
##                        or reference-skip spanning it).
##   - `footprint_depth`: count of reads that span this position at all,
##                        deletion/skip or not -- i.e. samtools mpileup's own
##                        per-site depth column.
## depth == footprint_depth means every spanning read has a matched base
## here. Any position where depth < footprint_depth (most notably depth == 0
## while footprint_depth > 0) means reads span it only via a deletion/skip;
## footprint_depth == 0 means no read spans the position at all. Left to
## downstream analysis to interpret as needed.
##
## Both numbers come from a single `samtools mpileup` pass: footprint_depth
## is mpileup's depth column as-is, and depth is derived by parsing the
## per-site pileup string and subtracting deletion/skip placeholders ("*",
## "<", ">") -- after first stripping "^X" read-start markers, since the
## character right after "^" encodes mapping quality and can itself render
## as "*"/"<"/">" by coincidence. This replaces running samtools depth twice
## (once default, once with -J), which re-scanned the same BAM region twice
## for no reason.
##
## Output is run-length encoded like a BedGraph: consecutive positions with
## identical (depth, footprint_depth) are collapsed into one interval, rather
## than emitting one row per position.
##
## If Intervals is supplied, scope is restricted to it twice over, since
## mpileup's own region handling isn't reliable enough to trust alone:
##   1. The BAM is filtered to only reads overlapping Intervals first
##      (samtools view -L), so mpileup never sees out-of-region coverage.
##      This matters because mpileup's docs warn that even a single -a can
##      fall back to padding whole, otherwise-untouched contigs if the BAM
##      has any coverage outside the requested regions.
##   2. mpileup itself is further restricted via -l Intervals.
##   3. The final collapsed table is intersected with Intervals again
##      (bedtools intersect), clipping any run that still extends past a
##      region boundary and dropping anything that slipped through.
## Without Intervals, coverage is reported genome-wide (note: the per-base
## pileup then scales with reference length, not BAM size -- size DiskGB
## accordingly).

workflow AssemblyCoverage {
    input {
        File Bam
        File BamIndex
        File? Intervals
        String Prefix
        String ImageTag = "latest"
    }

    call GetCoverage {
        input:
            Bam = Bam,
            BamIndex = BamIndex,
            Intervals = Intervals,
            Prefix = Prefix,
            ImageTag = ImageTag
    }

    output {
        File CoverageTable = GetCoverage.CoverageTable
    }
}

task GetCoverage {
    input {
        File Bam
        File BamIndex
        File? Intervals
        String Prefix
        String ImageTag = "latest"
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(Bam, "GB") * 2) + 20

    command <<<
        set -euo pipefail

        ln -s ~{Bam} input.bam
        ln -s ~{BamIndex} input.bam.bai

        # Collapses per-site mpileup rows into a BedGraph-like interval table. Defined
        # once and reused so the Intervals-filtering logic below doesn't have to
        # duplicate it; takes no File inputs itself, so it carries no localization
        # concerns of its own.
        collapse_pileup() {
            awk -F'\t' '
            BEGIN {
                OFS = "\t"
                print "chrom", "start", "end", "depth", "footprint_depth"
                first = 1
            }
            {
                cur_chrom = $1; pos = $2; footprint = $4; bases = $5
                gsub(/\^./, "", bases)
                placeholders = gsub(/[*<>]/, "", bases)
                d = footprint - placeholders

                if (!first && cur_chrom == chrom && d == depth && footprint == fdepth && pos == end + 1) {
                    end = pos
                } else {
                    if (!first) print chrom, start - 1, end, depth, fdepth
                    chrom = cur_chrom; start = pos; end = pos; depth = d; fdepth = footprint
                    first = 0
                }
            }
            END {
                if (!first) print chrom, start - 1, end, depth, fdepth
            }'
        }

        # Intervals is referenced only via the "+" placeholder idiom (never assigned to a
        # separate String/File declaration outside command) so Cromwell substitutes the
        # localized path here; converting it to a plain String earlier would keep the raw
        # gs:// URI and samtools would fail to read it.
        #
        # -B: skip BAQ recalibration (irrelevant without -f, and we don't use quality).
        # -Q 0 -q 0 -d 0: no base-quality, mapping-quality, or depth-cap filtering, so
        # this matches the unfiltered counting samtools depth did before.
        if ~{defined(Intervals)}; then
            samtools view -b ~{"-L " + Intervals} input.bam > subset.bam
            samtools index subset.bam

            samtools mpileup -a -B -Q 0 -q 0 -d 0 ~{"-l " + Intervals} subset.bam \
                | collapse_pileup > raw_coverage.tsv

            { head -n 1 raw_coverage.tsv
              tail -n +2 raw_coverage.tsv | bedtools intersect -a - ~{"-b " + Intervals}
            } > ~{Prefix}.coverage.tsv
        else
            samtools mpileup -a -a -B -Q 0 -q 0 -d 0 input.bam \
                | collapse_pileup > ~{Prefix}.coverage.tsv
        fi
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-bioinformatics:" + ImageTag
        memory: "8G"
        cpu: 2
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File CoverageTable = "~{Prefix}.coverage.tsv"
    }
}
