version 1.0

## Combines a reference-SVA bed and an MEI bed into one variant bed, then adds
## symmetric flanking on each side (clamped to chromosome bounds via
## `bedtools slop -g`), while preserving each variant's original (pre-flank)
## coordinates and strand as extra trailing columns for downstream distance
## calculations (e.g. scripts/process_methylation_window.R).
##
## Both input beds are expected to carry the variant ID in column 4 and
## strand in column 6 (SVA: chrom, start, end, sva_id, score, strand,
## subfamily; MEI: chrom, start, end, var_id, score, strand, min_length,
## max_length) -- everything else (score, subfamily, min/max length) is
## dropped, since only id/strand/coordinates are needed past this point.
## `#`-prefixed header lines in either input are skipped.

workflow CombineAndWindowVariants {
    input {
        File SVABed
        File MEIBed
        File ChromSizes
        Int FlankSize = 10000
        String Prefix = "combined"
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
    }

    call CombineAndWindow {
        input:
            SVABed = SVABed,
            MEIBed = MEIBed,
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
        File ChromSizes
        Int FlankSize = 10000
        String Prefix = "combined"
        String DockerImage
        Int MemoryGB = 4
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([SVABed, MEIBed], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        # Normalize both beds to a common 7-column layout -- chrom, start, end,
        # var_id, strand, var_start, var_end -- by column position (id@4,
        # strand@6 in both input schemas), duplicating start/end into
        # var_start/var_end *before* flanking so `bedtools slop` below only
        # ever touches columns 2/3, leaving the pre-flank coordinates in
        # columns 6/7 untouched.
        awk -F'\t' -v OFS='\t' '!/^#/ { print $1, $2, $3, $4, $6, $2, $3 }' ~{SVABed} ~{MEIBed} \
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
