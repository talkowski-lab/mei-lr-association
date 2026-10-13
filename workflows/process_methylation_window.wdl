version 1.0

## Bins per-CpG methylation calls from a BedtoolsIntersect output (ONE
## pb-CpG-tools pileup bed -- haplotype-specific or combined/unphased,
## doesn't matter which -- intersected against a window-expanded variant
## BED) into fixed 100bp-wide tiles by distance from the variant, and
## summarizes site count, coverage, and methylation score per (variant,
## haplotype, distance bin), via scripts/process_methylation_window.R.
## Haplotype is read from the pbcpg bed's own embedded column, not assumed
## from which file was passed in -- callers that need a per-haplotype
## breakdown call this once per haplotype file (see
## sva_mei_methylation_windows.wdl) and concatenate the resulting summaries
## afterward.

workflow ProcessMethylationWindow {
    input {
        String Individual
        File IntersectBed
        Int BinSize = 100
        Int FlankSize = 10000
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
        Int BinSize = 100
        Int FlankSize = 10000
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
            --indiv ~{Individual} \
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
