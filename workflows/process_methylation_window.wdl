version 1.0

## Bins per-5mC-site methylation calls from a BedtoolsIntersect output
## (5mC bedgraph intersected against a window-expanded variant BED) by
## distance from the variant, and summarizes site count and score
## statistics per (variant, distance bin), via
## scripts/process_methylation_window.R.

workflow ProcessMethylationWindow {
    input {
        String Individual
        File IntersectBed
        String ImageTag = "latest"
    }

    call SummarizeMethylationWindow {
        input:
            Individual = Individual,
            IntersectBed = IntersectBed,
            ImageTag = ImageTag
    }

    output {
        File MethylationWindowSummary = SummarizeMethylationWindow.MethylationWindowSummary
    }
}

task SummarizeMethylationWindow {
    input {
        String Individual
        File IntersectBed
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_gb = ceil(size(IntersectBed, "GB")) * 2 + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/process_methylation_window.R \
            --intersect-bed ~{IntersectBed} \
            --indiv ~{Individual} \
            --prefix ~{Individual}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: 2
        disks: "local-disk " + select_first([DiskGB, auto_disk_gb]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File MethylationWindowSummary = "~{Individual}.methylation_window_summary.tsv"
    }
}
