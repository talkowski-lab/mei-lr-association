version 1.0

## Combines non-reference (nrMEI) and reference (rMEI) SVA repeat-length data
## into a tidy per-haplotype length table plus per-individual length and
## deviance-from-mean matrices (raw and longest-pure-segment/LPS-adjusted).
## NrMEISvaRepeatSummary (ParseLongestPureSegment output for the nrMEI
## sequences) supplies the nrMEI LPS lengths; together with the optional
## RMEISvaRepeatSummary (ProcessRefSVAData's rMEISVARepeatSummaryFiltered) it
## is also written as a tidy per-haplotype repeat summary table. Driven by
## scripts/process_MEI_sva_data.R; nrMEIGtTable is expected to be the
## clustered GT table produced by cluster_mei_gt.R.

workflow ProcessMEISVAData {
    input {
        File NrMEIGtTable
        File NrMEISvaLength
        File NrMEISvaRepeatSummary
        File RMEISvaLengthIndiv
        File RMEISvaLengthLpsIndiv
        File RefSVALengthTable
        File RefSVABed
        File? RMEISvaRepeatSummary
        String ImageTag = "latest"
    }

    call ProcessMEISVADataTask {
        input:
            NrMEIGtTable = NrMEIGtTable,
            NrMEISvaLength = NrMEISvaLength,
            NrMEISvaRepeatSummary = NrMEISvaRepeatSummary,
            RMEISvaLengthIndiv = RMEISvaLengthIndiv,
            RMEISvaLengthLpsIndiv = RMEISvaLengthLpsIndiv,
            RefSVALengthTable = RefSVALengthTable,
            RefSVABed = RefSVABed,
            RMEISvaRepeatSummary = RMEISvaRepeatSummary,
            ImageTag = ImageTag
    }

    output {
        File MEILengthsTidy = ProcessMEISVADataTask.MEILengthsTidy
        File MEILengthPerIndivMatrix = ProcessMEISVADataTask.MEILengthPerIndivMatrix
        File MEIDeviancePerIndivMatrix = ProcessMEISVADataTask.MEIDeviancePerIndivMatrix
        File MEILengthLpsPerIndivMatrix = ProcessMEISVADataTask.MEILengthLpsPerIndivMatrix
        File MEIDevianceLpsPerIndivMatrix = ProcessMEISVADataTask.MEIDevianceLpsPerIndivMatrix
        File MEIRepeatSummaryTidy = ProcessMEISVADataTask.MEIRepeatSummaryTidy
    }
}

task ProcessMEISVADataTask {
    input {
        File NrMEIGtTable
        File NrMEISvaLength
        File NrMEISvaRepeatSummary
        File RMEISvaLengthIndiv
        File RMEISvaLengthLpsIndiv
        File RefSVALengthTable
        File RefSVABed
        File? RMEISvaRepeatSummary
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil((size([NrMEIGtTable, NrMEISvaLength, NrMEISvaRepeatSummary, RMEISvaLengthIndiv, RMEISvaLengthLpsIndiv, RefSVALengthTable, RefSVABed], "GB") + size(RMEISvaRepeatSummary, "GB")) * 2) + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/process_MEI_sva_data.R \
            --nrMEI-gt-table ~{NrMEIGtTable} \
            --nrMEI-sva-length ~{NrMEISvaLength} \
            --nrMEI-sva-repeat-summary ~{NrMEISvaRepeatSummary} \
            --rMEI-sva-length-indiv ~{RMEISvaLengthIndiv} \
            --rMEI-sva-length-lps-indiv ~{RMEISvaLengthLpsIndiv} \
            --refSVA-length-table ~{RefSVALengthTable} \
            --refSVA-bed ~{RefSVABed} \
            ~{"--rMEI-sva-repeat-summary " + RMEISvaRepeatSummary}
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
        File MEILengthsTidy = "SVA_MEI_lengths_tidy.txt.gz"
        File MEILengthPerIndivMatrix = "mei_SVA_length_perindiv_matrix.bed.gz"
        File MEIDeviancePerIndivMatrix = "mei_SVA_deviance_perindiv_matrix.bed.gz"
        File MEILengthLpsPerIndivMatrix = "mei_SVA_length_lps_perindiv_matrix.bed.gz"
        File MEIDevianceLpsPerIndivMatrix = "mei_SVA_deviance_lps_perindiv_matrix.bed.gz"
        File MEIRepeatSummaryTidy = "SVA_MEI_repeat_summary_tidy.txt.gz"
    }
}
