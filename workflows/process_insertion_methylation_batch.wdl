version 1.0

import "utils/concat_files.wdl" as ConcatFiles

## Batch variant of process_insertion_methylation.wdl: summarizes per-read
## insertion methylation calls against MEI genotypes for an array of
## individuals and their corresponding methylation tables.
##
## Individuals/MethylationTables are chunked into batches of BatchSize (all
## of them, in one batch, if BatchSize is unset). Each batch is handed to
## its own scatter shard, which runs SummarizeInsertionMethylationBatch --
## looping over that batch's individuals sequentially on a single VM,
## rather than via `scatter` per individual, since each individual's run is
## short enough that per-individual VM spin-up would dominate the cost --
## and then concatenates that shard's per-individual results before the
## shard finishes. The shards' results are then concatenated together into
## the final combined outputs.
##
## Each individual's rows are tagged with `indiv` inside
## process_insertion_methylation.R itself, so the per-batch task can gather
## its per-individual output files via glob() (the file-array pattern this
## repo's backend reliably delocalizes, e.g. render_rmarkdown.wdl) without
## needing glob's returned order to match Individuals order.
##
## GenotypeTable and MeiReferenceTable are already combined across
## SVA/LINE1/Alu and cover all individuals; they're shared, constant inputs
## across every batch.

workflow ProcessInsertionMethylationBatch {
    input {
        Array[String] Individuals
        Array[File] MethylationTables
        File GenotypeTable
        File MeiReferenceTable
        String Prefix
        Int? BatchSize
        String ImageTag = "latest"
        Int MinInsertionLength = 20
        Int MatchBp = 10
        Int MatchPerc = 3
    }

    Int total_individuals = length(Individuals)
    Int effective_batch_size = select_first([BatchSize, total_individuals])
    Int n_batches = (total_individuals + effective_batch_size - 1) / effective_batch_size

    scatter (batch_idx in range(n_batches)) {
        # No slicing operator in WDL 1.0, so Individuals[batch_idx*size : ...]
        # isn't available -- gather this batch's indices via select_all over a
        # masked range, then gather the actual elements by those indices.
        scatter (pos in range(total_individuals)) {
            if (pos >= batch_idx * effective_batch_size && pos < (batch_idx + 1) * effective_batch_size) {
                Int keep_pos = pos
            }
        }
        Array[Int] batch_indices = select_all(keep_pos)

        scatter (idx in batch_indices) {
            String batch_indiv = Individuals[idx]
            File batch_methyl_table = MethylationTables[idx]
        }

        call SummarizeInsertionMethylationBatch {
            input:
                Individuals = batch_indiv,
                MethylationTables = batch_methyl_table,
                GenotypeTable = GenotypeTable,
                MeiReferenceTable = MeiReferenceTable,
                ImageTag = ImageTag,
                MinInsertionLength = MinInsertionLength,
                MatchBp = MatchBp,
                MatchPerc = MatchPerc
        }

        call ConcatFiles.ConcatenateDelim as ConcatBatchDataSumm {
            input:
                InputFiles = SummarizeInsertionMethylationBatch.DataSummFiles,
                OutputName = "batch_" + batch_idx + ".ins_methyl_data_summ.tsv",
                ImageTag = ImageTag
        }

        call ConcatFiles.ConcatenateDelim as ConcatBatchHapSumm {
            input:
                InputFiles = SummarizeInsertionMethylationBatch.HapSummFiles,
                OutputName = "batch_" + batch_idx + ".ins_methyl_hap_summ.tsv",
                ImageTag = ImageTag
        }
    }

    # Batch outputs already have (at most) one header each and their rows are
    # already indiv-tagged (by the R script), so this call only
    # re-concatenates across batches.
    call ConcatFiles.ConcatenateDelim as MergeDataSumm {
        input:
            InputFiles = ConcatBatchDataSumm.ConcatenatedFile,
            OutputName = Prefix + ".ins_methyl_data_summ.tsv",
            ImageTag = ImageTag
    }

    call ConcatFiles.ConcatenateDelim as MergeHapSumm {
        input:
            InputFiles = ConcatBatchHapSumm.ConcatenatedFile,
            OutputName = Prefix + ".ins_methyl_hap_summ.tsv",
            ImageTag = ImageTag
    }

    output {
        File InsMethylDataSumm = MergeDataSumm.ConcatenatedFile
        File InsMethylHapSumm = MergeHapSumm.ConcatenatedFile
    }
}

task SummarizeInsertionMethylationBatch {
    input {
        Array[String] Individuals
        Array[File] MethylationTables
        File GenotypeTable
        File MeiReferenceTable
        String ImageTag = "latest"
        Int MinInsertionLength = 20
        Int MatchBp = 10
        Int MatchPerc = 3
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_gb = ceil(size(MethylationTables, "GB") + size(GenotypeTable, "GB") + size(MeiReferenceTable, "GB")) * 2 + 10

    command <<<
        set -euo pipefail

        mkdir -p outputs

        indivs=(~{sep=" " Individuals})
        methyl_tables=(~{sep=" " MethylationTables})

        # Sequential, single-VM loop rather than a `scatter` -- each
        # individual's run is short enough that per-individual VM spin-up
        # would dominate the cost. Each individual's rows are tagged with
        # its `indiv` inside the R script itself, so the outputs below can
        # be gathered by glob() without depending on file order matching
        # Individuals order.
        for i in "${!indivs[@]}"; do
            indiv="${indivs[$i]}"
            methyl_table="${methyl_tables[$i]}"

            Rscript /scripts/process_insertion_methylation.R \
                --indiv "$indiv" \
                --methylation-table "$methyl_table" \
                --genotype-table ~{GenotypeTable} \
                --mei-reference-table ~{MeiReferenceTable} \
                --min-insertion-length ~{MinInsertionLength} \
                --match-bp ~{MatchBp} \
                --match-perc ~{MatchPerc} \
                --prefix "outputs/$indiv"
        done
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
        Array[File] DataSummFiles = glob("outputs/*.ins_methyl_data_summ.tsv")
        Array[File] HapSummFiles = glob("outputs/*.ins_methyl_hap_summ.tsv")
    }
}
