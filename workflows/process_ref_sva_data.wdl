version 1.0

## Combines per-haplotype ref-SVA repeat-length tables with QC flags and a
## cohort SV VCF to produce per-individual repeat-length matrices (filtered
## to well-assembled SVAs), a per-haplotype rMEI genotype table in the nrMEI
## tidy layout, plus a genotype-vs-assembly concordance check
## against MED (mobile element deletion) calls from L1-MEAID. Reference SVAs
## matched to a MED call are polymorphic, the rest fixed; the filtered
## per-haplotype tables are split into rMEI_SVA_* and fixed_SVA_*, the
## per-individual matrices cover both. RefSVARepeatSummaryTable (the
## ParseLongestPureSegment output, with ID already split into indiv/hap/sva_id)
## supplies the longest-pure-segment lengths and the per-haplotype repeat
## summary table. Driven by
## scripts/process_ref_sva_data.R.

workflow ProcessRefSVAData {
    input {
        File RefSVALengthTable
        File RefSVAFlagTable
        File RefSVABed
        String SVVCFUrl
        File L1MEAIDTable
        File RefSVARepeatSummaryTable
        String ImageTag = "latest"
    }

    call ProcessRefSVADataTask {
        input:
            RefSVALengthTable = RefSVALengthTable,
            RefSVAFlagTable = RefSVAFlagTable,
            RefSVABed = RefSVABed,
            SVVCFUrl = SVVCFUrl,
            L1MEAIDTable = L1MEAIDTable,
            RefSVARepeatSummaryTable = RefSVARepeatSummaryTable,
            ImageTag = ImageTag
    }

    output {
        File FixedSVALengthFiltered = ProcessRefSVADataTask.FixedSVALengthFiltered
        File rMEISVALengthFiltered = ProcessRefSVADataTask.rMEISVALengthFiltered
        File rMEISVATotLengthPerIndiv = ProcessRefSVADataTask.rMEISVATotLengthPerIndiv
        File rMEISVATotLengthLpsPerIndiv = ProcessRefSVADataTask.rMEISVATotLengthLpsPerIndiv
        File RefSVALengthPerIndivMatrix = ProcessRefSVADataTask.RefSVALengthPerIndivMatrix
        File RefSVALengthLpsPerIndivMatrix = ProcessRefSVADataTask.RefSVALengthLpsPerIndivMatrix
        File FixedSVARepeatSummaryFiltered = ProcessRefSVADataTask.FixedSVARepeatSummaryFiltered
        File rMEISVARepeatSummaryFiltered = ProcessRefSVADataTask.rMEISVARepeatSummaryFiltered
        File rMEISVAGtTidy = ProcessRefSVADataTask.rMEISVAGtTidy
    }
}

task ProcessRefSVADataTask {
    input {
        File RefSVALengthTable
        File RefSVAFlagTable
        File RefSVABed
        String SVVCFUrl
        File L1MEAIDTable
        File RefSVARepeatSummaryTable
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([RefSVALengthTable, RefSVARepeatSummaryTable, RefSVAFlagTable, RefSVABed, L1MEAIDTable], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/process_ref_sva_data.R \
            --ref-sva-length-table ~{RefSVALengthTable} \
            --ref-sva-flag-table ~{RefSVAFlagTable} \
            --ref-sva-bed ~{RefSVABed} \
            --sv-vcf ~{SVVCFUrl} \
            --l1meaid-table ~{L1MEAIDTable} \
            --ref-sva-repeat-summary-table ~{RefSVARepeatSummaryTable}
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
        File FixedSVALengthFiltered = "fixed_SVA_length_filtered.txt.gz"
        File rMEISVALengthFiltered = "rMEI_SVA_length_filtered.txt.gz"
        File rMEISVATotLengthPerIndiv = "rMEI_SVA_totlength_perindiv_filtered.txt.gz"
        File rMEISVATotLengthLpsPerIndiv = "rMEI_SVA_totlength_lps_perindiv_filtered.txt.gz"
        File RefSVALengthPerIndivMatrix = "ref_SVA_length_perindiv_matrix.bed.gz"
        File RefSVALengthLpsPerIndivMatrix = "ref_SVA_length_lps_perindiv_matrix.bed.gz"
        File FixedSVARepeatSummaryFiltered = "fixed_SVA_repeat_summary_filtered.txt.gz"
        File rMEISVARepeatSummaryFiltered = "rMEI_SVA_repeat_summary_filtered.txt.gz"
        File rMEISVAGtTidy = "rMEI_SVA_gt_tidy.txt.gz"
    }
}
