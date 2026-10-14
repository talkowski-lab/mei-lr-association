version 1.0

## Bins per-CpG methylation calls from a BedtoolsIntersect output (one
## pb-CpG-tools pileup bed intersected against a window-expanded variant
## BED) into fixed-width tiles by distance from the variant, and summarizes
## per (variant, haplotype, distance bin), via
## scripts/process_methylation_window.R.

workflow ProcessMethylationWindow {
    input {
        String Individual
        File IntersectBed
        Int BinSize = 200
        Int FlankSize = 5000
        String Prefix = Individual
        String ImageTag = "latest"
    }

    call SummarizeMethylationWindow {
        input:
            Individual = Individual,
            IntersectBed = IntersectBed,
            BinSize = BinSize,
            FlankSize = FlankSize,
            Prefix = Prefix,
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
        Int BinSize = 200
        Int FlankSize = 5000
        String Prefix = Individual
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_gb = ceil(size(IntersectBed, "GB")) * 2 + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/process_methylation_window.R \
            --intersect-bed ~{IntersectBed} \
            --bin-size ~{BinSize} \
            --flank-size ~{FlankSize} \
            --prefix ~{Prefix}
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
        File MethylationWindowSummary = "~{Prefix}.methylation_window_summary.tsv"
    }
}
