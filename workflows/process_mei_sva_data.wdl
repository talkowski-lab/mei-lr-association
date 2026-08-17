version 1.0

## Combines non-reference (nrMEI) and reference (rMEI) SVA repeat-length data
## into a tidy per-haplotype length table plus per-individual length and
## deviance-from-mean matrices (raw and longest-pure-segment/LPS-adjusted).
## Driven by scripts/process_MEI_sva_data.R; nrMEIGtTable is expected to be
## the clustered GT table produced by cluster_mei_gt.R.

workflow ProcessMEISVAData {
    input {
        File NrMEIGtTable
        File NrMEISvaLength
        File NrMEISvaLengthLps
        File RMEISvaLengthIndiv
        File RMEISvaLengthLpsIndiv
        File RefSVALengthTable
        File RefSVABed
        String ImageTag = "latest"
    }

    call ProcessMEISVADataTask {
        input:
            NrMEIGtTable = NrMEIGtTable,
            NrMEISvaLength = NrMEISvaLength,
            NrMEISvaLengthLps = NrMEISvaLengthLps,
            RMEISvaLengthIndiv = RMEISvaLengthIndiv,
            RMEISvaLengthLpsIndiv = RMEISvaLengthLpsIndiv,
            RefSVALengthTable = RefSVALengthTable,
            RefSVABed = RefSVABed,
            ImageTag = ImageTag
    }

    output {
        File MEILengthsTidy = ProcessMEISVADataTask.MEILengthsTidy
        File MEILengthPerIndivMatrix = ProcessMEISVADataTask.MEILengthPerIndivMatrix
        File MEIDeviancePerIndivMatrix = ProcessMEISVADataTask.MEIDeviancePerIndivMatrix
        File MEILengthLpsPerIndivMatrix = ProcessMEISVADataTask.MEILengthLpsPerIndivMatrix
        File MEIDevianceLpsPerIndivMatrix = ProcessMEISVADataTask.MEIDevianceLpsPerIndivMatrix
    }
}

task ProcessMEISVADataTask {
    input {
        File NrMEIGtTable
        File NrMEISvaLength
        File NrMEISvaLengthLps
        File RMEISvaLengthIndiv
        File RMEISvaLengthLpsIndiv
        File RefSVALengthTable
        File RefSVABed
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([NrMEIGtTable, NrMEISvaLength, NrMEISvaLengthLps, RMEISvaLengthIndiv, RMEISvaLengthLpsIndiv, RefSVALengthTable, RefSVABed], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/process_MEI_sva_data.R \
            --nrMEI-gt-table ~{NrMEIGtTable} \
            --nrMEI-sva-length ~{NrMEISvaLength} \
            --nrMEI-sva-length-lps ~{NrMEISvaLengthLps} \
            --rMEI-sva-length-indiv ~{RMEISvaLengthIndiv} \
            --rMEI-sva-length-lps-indiv ~{RMEISvaLengthLpsIndiv} \
            --refSVA-length-table ~{RefSVALengthTable} \
            --refSVA-bed ~{RefSVABed} 
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
    }
}
