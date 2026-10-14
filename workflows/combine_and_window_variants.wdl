version 1.0

## Combines a reference-SVA bed, an MEI bed, and a third bed (same format as
## the MEI bed) into one variant bed, dedupes exact-duplicate rows, and
## flanks by FlankSize on each side via `bedtools slop -g`. Each bed's first
## 6 columns (chrom, start, end, id, score, strand) are taken as-is; anything
## past column 6 is dropped. Output: chrom, start, end, id, score, strand,
## var_start, var_end (pre-flank start/end).

workflow CombineAndWindowVariants {
    input {
        File SVABed
        File MEIBed
        File OtherBed
        File ChromSizes
        Int FlankSize = 5000
        String Prefix = "combined"
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
    }

    call CombineAndWindow {
        input:
            SVABed = SVABed,
            MEIBed = MEIBed,
            OtherBed = OtherBed,
            ChromSizes = ChromSizes,
            FlankSize = FlankSize,
            Prefix = Prefix,
            DockerImage = DockerImage
    }

    output {
        File WindowedVariantBed = CombineAndWindow.WindowedVariantBed
    }
}

task CombineAndWindow {
    input {
        File SVABed
        File MEIBed
        File OtherBed
        File ChromSizes
        Int FlankSize = 5000
        String Prefix = "combined"
        String DockerImage
        Int MemoryGB = 4
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([SVABed, MEIBed, OtherBed], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        awk -F'\t' -v OFS='\t' '
            !/^#/ {
                if (NF < 6) {
                    printf "%s line %d: expected >=6 tab-separated columns (chrom, start, end, id, score, strand), got %d: %s\n", FILENAME, FNR, NF, $0 > "/dev/stderr"
                    exit 1
                }
                print $1, $2, $3, $4, $5, $6, $2, $3
            }
        ' ~{SVABed} ~{MEIBed} ~{OtherBed} \
            | sort -u \
            | sort -k1,1 -k2,2n > combined.bed

        bedtools slop -b ~{FlankSize} -i combined.bed -g ~{ChromSizes} > "~{Prefix}.windowed_variants.bed"
    >>>

    runtime {
        docker: DockerImage
        memory: MemoryGB + " GB"
        cpu: 1
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File WindowedVariantBed = "~{Prefix}.windowed_variants.bed"
    }
}
