version 1.0

# Filter the Assoc Table to those with a pvalue below PValueThreshold, then
# join each significant (sva_id, feature) pair against the SVA genotype/
# repeat-length BED to produce the list of significant pairs plus one wide
# row per (variant, feature) combination. Also drops any feature not present
# in the first column of AVTable (it would otherwise fail a per-feature
# AVTable lookup downstream), reporting dropped features to stderr.
# If SignatureType is splice, the leafcutter feature is collapsed to its gene id.
task SubsetSignificantSVAGT {
  input {
    String SignatureType
    File QTLAssocSummary
    File SVAGTBed
    File AVTable
    Float PValueThreshold = 0.01
    String Prefix = "sig_sva"
    String ImageTag = "latest"
    Int MemoryGB = 4
    Int? DiskGB
  }

  Int auto_disk_size = ceil(size([QTLAssocSummary, SVAGTBed, AVTable], "GB") * 2) + 10

  command <<<
    set -euo pipefail

    Rscript /scripts/subset_significant_sva_gt.R \
        --qtl-assoc-summary ~{QTLAssocSummary} \
        --sva-gt-bed ~{SVAGTBed} \
        --avtable-file ~{AVTable} \
        --signature-type ~{SignatureType} \
        --p-threshold ~{PValueThreshold} \
        --prefix ~{Prefix}
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
    File SignificantPairs = "~{Prefix}_sig_pairs.txt"
    File SignificantSVAGT = "~{Prefix}_sig_sva_gt.txt.gz"
    Array[String] FeatureNames = read_json("feature_names.json")
    Array[Array[String]] SVAIDs = read_json("sva_ids.json")
  }
}

# Looks up one row (by GeneName, matched against the table's first/id
# column) in an already-downloaded QTL data table and pulls the requested
# fine-mapping input columns.
task GatherSNVFM_Input {

  input {
    File AVTableFile
    String GeneName
    String ImageTag = "latest"
    Int MemoryGB = 4
    Int? DiskGB
  }

  command <<<
    set -euo pipefail

    Rscript /scripts/gather_snv_fm_input.R \
        --avtable-file ~{AVTableFile} \
        --gene-name ~{GeneName}
  >>>

  runtime {
    docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
    memory: MemoryGB + " GB"
    cpu: 2
    disks: "local-disk " + select_first([DiskGB, 20]) + " SSD"
    preemptible: 3
    maxRetries: 2
  }

  output {
    String GenotypeDosages = read_string("genotype_dosages.txt")
    String GenotypeDosagesIndex = read_string("genotype_dosages_index.txt")
    String TensorQTLPermutations = read_string("tensorqtl_permutations.txt")
    String PhenotypeBed = read_string("phenotype_bed.txt")
  }
}

# This is a straight up copy of
# https://raw.githubusercontent.com/AoU-Multiomics-Analysis/susieR/refs/heads/main/workflows/susieRonly.wdl
# On 8/19/2026
# Except for the hard disk requirement because 500Gb is too much

task susieR {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        File PhenotypeBed
        Int CisDistance
        String OutputPrefix
        Int memory
        Int NumPrempt
        Float? MAF
        Boolean MatchPhenotypeIDSubstring = false
        Boolean ReuseGenotypeMatrix = false
        Boolean SelectTopPhenotypePerCluster = false
        String TopPhenotypePerClusterPvalueColumn = "qval"
        Int PreparedWindowSize = -1
        File? VariantList
        File? AncestryFile
        File? AdditionalGenotypesBed
        Int? DiskGB
    }
    String phenotype_match_mode = if MatchPhenotypeIDSubstring then "contains" else "exact"

    # Genotype dosages tend to dominate input size; susie's outputs (full
    # posterior parquet, RDS object, etc.) can be a comparable multiple of
    # that, so weight by 3x rather than the flat 2x used elsewhere in this
    # repo, plus a 20 GB floor for small inputs.
    Int auto_disk_size = ceil(
        size([GenotypeDosages, GenotypeDosageIndex, QTLCovariates, TensorQTLPermutations, SampleList, PhenotypeBed], "GB") * 3
        + size(AdditionalGenotypesBed, "GB") * 3
        + size(VariantList, "GB")
        + size(AncestryFile, "GB")
    ) + 20

    command <<<

        echo "Phenotype match mode: ~{phenotype_match_mode}"
        echo "Output prefix selected for this run: ~{OutputPrefix}"
        zcat ~{PhenotypeBed} | head -n 1 > header.txt
        if [ "~{phenotype_match_mode}" = "contains" ]; then
            zcat ~{PhenotypeBed} \
                | awk -v needle="~{OutputPrefix}" 'FNR == 1 && $4 == "phenotype_id" {next} index($4, needle) > 0' \
                > input_gene.txt
        else
            zcat ~{PhenotypeBed} \
                | awk -v phenotype_id="~{OutputPrefix}" 'FNR == 1 && $4 == "phenotype_id" {next} $4 == phenotype_id' \
                > input_gene.txt
        fi
        if [ ! -s input_gene.txt ]; then
            echo "No rows in PhenotypeBed matched the requested phenotype IDs" >&2
            exit 1
        fi
        head -n 1 input_gene.txt | awk -F'\t' 'BEGIN{OFS="\t"} {$4="skip"; print}' > skip.txt        
        cat header.txt input_gene.txt skip.txt > input_gene.bed  

        if [ ~{PreparedWindowSize} -ge 0 ] && [ ~{CisDistance} -gt ~{PreparedWindowSize} ]; then
            echo "CisDistance (~{CisDistance}) is larger than the prepared dosage WindowSize (~{PreparedWindowSize}); prepared dosages may not contain all requested variants." >&2
            echo "Re-run PrepInputs with WindowSize >= CisDistance or reduce CisDistance." >&2
            exit 1
        fi

        Rscript /tmp/susie.R ~{if defined(MAF) then "--MAF ~{MAF}  " else ""} ~{if ReuseGenotypeMatrix then "--reuse_genotype_matrix true  " else ""} ~{if SelectTopPhenotypePerCluster then "--select_top_phenotype_per_cluster true --top_phenotype_pvalue_column " else ""}~{if SelectTopPhenotypePerCluster then TopPhenotypePerClusterPvalueColumn else ""} ~{if defined(AncestryFile) then "--AncestryMetadata ~{AncestryFile}  "  else ""} ~{if defined(VariantList) then "--VariantList ~{VariantList}  "  else ""}  ~{if defined(AdditionalGenotypesBed) then "--AdditionalGenotypesBed ~{AdditionalGenotypesBed}  "  else ""} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix input_gene.bed \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance} 

    >>>

    runtime {
        docker: 'ghcr.io/aou-multiomics-analysis/susier:main'
        memory: "${memory}GB"
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        bootDiskSizeGb: 25
        preemptible: "${NumPrempt}"
        cpu: 2
    }

    output {
        File SusieParquet = "${OutputPrefix}.parquet"
        File lbfParquet = "${OutputPrefix}.lbf_variable.parquet"
        File FullSusieParquet = "${OutputPrefix}.full_susie.parquet"
        File VariantPositionSummary = "${OutputPrefix}.variant_position_summary.tsv"
        File SusieObject = "${OutputPrefix}_susie.rds"
    }
}
