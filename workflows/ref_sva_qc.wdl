version 1.0

## Annotates reference SVA elements with QC flags, combining:
##   - N-content in the extracted ref SVA sequences (per haplotype)
##   - assembly coverage status over each SVA's footprint (per haplotype),
##     flagging elements with zero coverage ("Missing") or multi-mapping
##     coverage ("Remap"/"Misc")
##   - large skips (deletions relative to the repeat consensus) within the
##     VNTR/hexamer repeat regions of each SVA
##
## Produces a flag table (one row per SVA x haplotype) and a BED of SVA
## regions worth re-mapping, driven by scripts/ref_sva_qc.R.

workflow RefSvaQc {
    input {
        File SvaRefBed
        File SvaRefSeqs
        File SvaRefCovTableHap1
        File SvaRefCovTableHap2
        File SvaRefRepCoords
        String Prefix
        String ImageTag = "latest"
    }

    call AnnotateRefSvaFlags {
        input:
            SvaRefBed = SvaRefBed,
            SvaRefSeqs = SvaRefSeqs,
            SvaRefCovTableHap1 = SvaRefCovTableHap1,
            SvaRefCovTableHap2 = SvaRefCovTableHap2,
            SvaRefRepCoords = SvaRefRepCoords,
            Prefix = Prefix,
            ImageTag = ImageTag
    }

    output {
        File RemapRegionsBed = AnnotateRefSvaFlags.RemapRegionsBed
        File RefSvaFlags = AnnotateRefSvaFlags.RefSvaFlags
    }
}

task AnnotateRefSvaFlags {
    input {
        File SvaRefBed
        File SvaRefSeqs
        File SvaRefCovTableHap1
        File SvaRefCovTableHap2
        File SvaRefRepCoords
        String Prefix
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([SvaRefBed, SvaRefSeqs, SvaRefCovTableHap1, SvaRefCovTableHap2, SvaRefRepCoords], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/ref_sva_qc.R \
            --sva-ref-bed ~{SvaRefBed} \
            --sva-ref-seqs ~{SvaRefSeqs} \
            --sva-ref-cov-table-hap1 ~{SvaRefCovTableHap1} \
            --sva-ref-cov-table-hap2 ~{SvaRefCovTableHap2} \
            --sva-ref-rep-coords ~{SvaRefRepCoords} \
            --prefix ~{Prefix}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: 2
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File RemapRegionsBed = "~{Prefix}_remap_regions.bed"
        File RefSvaFlags = "~{Prefix}_ref_sva_flags.txt"
    }
}
