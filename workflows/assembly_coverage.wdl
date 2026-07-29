version 1.0

## Computes per-base coverage across an assembly-to-reference alignment BAM
## (e.g. contigs aligned to the reference with minimap2 asm5), reporting two
## depths side by side so gaps vs. deletions can be told apart downstream:
##   - `depth`:           default samtools depth behavior, which excludes
##                        deleted reference bases -- a read with a deletion
##                        contributes 0 here.
##   - `footprint_depth`: samtools depth -J, which counts a read even where
##                        it has a deletion, i.e. "does any read/contig span
##                        this position at all".
## depth == footprint_depth means the position is covered by matched bases
## the same way it's spanned. Any position where depth < footprint_depth
## (most notably depth == 0 while footprint_depth > 0) means reads span it
## but only via a deletion there; footprint_depth == 0 means no read spans
## the position at all. Left to downstream analysis to interpret as needed.
##
## Output is run-length encoded like a BedGraph: consecutive positions with
## identical (depth, footprint_depth) are collapsed into one interval, rather
## than emitting one row per position.
##
## If Intervals is supplied, both depth passes are restricted to those
## regions via samtools depth's own -b option; otherwise coverage is reported
## genome-wide (note: without Intervals the per-base output scales with
## reference length, not BAM size -- size DiskGB accordingly).

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

    Int auto_disk_size = ceil(size(Bam, "GB") * 4) + 20

    command <<<
        set -euo pipefail

        ln -s ~{Bam} input.bam
        ln -s ~{BamIndex} input.bam.bai

        # Intervals is referenced only via the "+" placeholder idiom (never assigned to a
        # separate String/File declaration outside command) so Cromwell substitutes the
        # localized path here; converting it to a plain String earlier would keep the raw
        # gs:// URI and samtools would fail to read it.
        samtools depth -a ~{"-b " + Intervals} input.bam > depth.tsv
        samtools depth -a -J ~{"-b " + Intervals} input.bam > footprint_depth.tsv

        # Collapse consecutive positions with identical (depth, footprint_depth) into a
        # single interval, like a BedGraph, instead of one row per position.
        paste depth.tsv footprint_depth.tsv | awk -F'\t' 'BEGIN {
            OFS = "\t"
            print "chrom", "start", "end", "depth", "footprint_depth"
            first = 1
        } {
            d = $3; f = $6
            if (!first && $1 == chrom && d == depth && f == footprint && $2 == end + 1) {
                end = $2
            } else {
                if (!first) print chrom, start - 1, end, depth, footprint
                chrom = $1; start = $2; end = $2; depth = d; footprint = f
                first = 0
            }
        }
        END {
            if (!first) print chrom, start - 1, end, depth, footprint
        }' > ~{Prefix}.coverage.tsv
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-coverage:" + ImageTag
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
