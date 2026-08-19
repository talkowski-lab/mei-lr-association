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
