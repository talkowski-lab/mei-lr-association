version 1.0

## Summarizes per-read insertion methylation calls (from
## ExtractInsertionMethylation) against MEI genotypes, for one individual
## and one methylation table, via scripts/process_insertion_methylation.R.
## GenotypeTable and MeiReferenceTable are already combined across
## SVA/LINE1/Alu and cover all individuals; this task filters internally to
## Individual. See process_insertion_methylation_batch.wdl for the
## single-VM, multi-individual variant.

workflow ProcessInsertionMethylation {
    input {
        String Individual
        File MethylationTable
        File GenotypeTable
        File MeiReferenceTable
        String ImageTag = "latest"
        Int MinInsertionLength = 20
        Int MatchBp = 10
        Int MatchPerc = 3
    }

    call SummarizeInsertionMethylation {
        input:
            Individual = Individual,
            MethylationTable = MethylationTable,
            GenotypeTable = GenotypeTable,
            MeiReferenceTable = MeiReferenceTable,
            ImageTag = ImageTag,
            MinInsertionLength = MinInsertionLength,
            MatchBp = MatchBp,
            MatchPerc = MatchPerc
    }

    output {
        File InsMethylDataSumm = SummarizeInsertionMethylation.InsMethylDataSumm
        File InsMethylHapSumm = SummarizeInsertionMethylation.InsMethylHapSumm
    }
}

task SummarizeInsertionMethylation {
    input {
        String Individual
        File MethylationTable
        File GenotypeTable
        File MeiReferenceTable
        String ImageTag = "latest"
        Int MinInsertionLength = 20
        Int MatchBp = 10
        Int MatchPerc = 3
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_gb = ceil(size(MethylationTable, "GB") + size(GenotypeTable, "GB") + size(MeiReferenceTable, "GB")) * 2 + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/process_insertion_methylation.R \
            --indiv ~{Individual} \
            --methylation-table ~{MethylationTable} \
            --genotype-table ~{GenotypeTable} \
            --mei-reference-table ~{MeiReferenceTable} \
            --min-insertion-length ~{MinInsertionLength} \
            --match-bp ~{MatchBp} \
            --match-perc ~{MatchPerc} \
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
        File InsMethylDataSumm = "~{Individual}.ins_methyl_data_summ.tsv"
        File InsMethylHapSumm = "~{Individual}.ins_methyl_hap_summ.tsv"
    }
}
