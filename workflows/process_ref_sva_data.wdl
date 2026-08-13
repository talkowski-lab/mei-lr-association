version 1.0

## Combines per-haplotype ref-SVA repeat-length tables with QC flags and a
## cohort SV VCF to produce per-individual repeat-length matrices (filtered
## to well-assembled SVAs), plus a genotype-vs-assembly concordance check
## against MED (mobile element deletion) calls from L1-MEAID. Driven by
## scripts/process_ref_sva_data.R.

workflow ProcessRefSVAData {
    input {
        File RefSVALengthTable
        File RefSVALpsTable
        File RefSVAFlagTable
        File RefSVABed
        String SVVCFUrl
        File L1MEAIDTable
        String ImageTag = "latest"
    }

    call ProcessRefSVADataTask {
        input:
            RefSVALengthTable = RefSVALengthTable,
            RefSVALpsTable = RefSVALpsTable,
            RefSVAFlagTable = RefSVAFlagTable,
            RefSVABed = RefSVABed,
            SVVCFUrl = SVVCFUrl,
            L1MEAIDTable = L1MEAIDTable,
            ImageTag = ImageTag
    }

    output {
        File RefSVALengthFiltered = ProcessRefSVADataTask.RefSVALengthFiltered
        File rMEISVATotLengthPerIndiv = ProcessRefSVADataTask.rMEISVATotLengthPerIndiv
        File rMEISVATotLengthLpsPerIndiv = ProcessRefSVADataTask.rMEISVATotLengthLpsPerIndiv
        File RefSVALengthPerIndivMatrix = ProcessRefSVADataTask.RefSVALengthPerIndivMatrix
        File RefSVALengthLpsPerIndivMatrix = ProcessRefSVADataTask.RefSVALengthLpsPerIndivMatrix
    }
}

task ProcessRefSVADataTask {
    input {
        File RefSVALengthTable
        File RefSVALpsTable
        File RefSVAFlagTable
        File RefSVABed
        String SVVCFUrl
        File L1MEAIDTable
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([RefSVALengthTable, RefSVALpsTable, RefSVAFlagTable, RefSVABed, L1MEAIDTable], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/process_ref_sva_data.R \
            --ref-sva-length-table ~{RefSVALengthTable} \
            --ref-sva-lps-table ~{RefSVALpsTable} \
            --ref-sva-flag-table ~{RefSVAFlagTable} \
            --ref-sva-bed ~{RefSVABed} \
            --sv-vcf ~{SVVCFUrl} \
            --l1meaid-table ~{L1MEAIDTable}
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
        File RefSVALengthFiltered = "ref_SVA_length_filtered.txt.gz"
        File rMEISVATotLengthPerIndiv = "rMEI_SVA_totlength_perindiv_filtered.txt.gz"
        File rMEISVATotLengthLpsPerIndiv = "rMEI_SVA_totlength_lps_perindiv_filtered.txt.gz"
        File RefSVALengthPerIndivMatrix = "ref_SVA_length_perindiv_matrix.bed.gz"
        File RefSVALengthLpsPerIndivMatrix = "ref_SVA_length_lps_perindiv_matrix.bed.gz"
    }
}
