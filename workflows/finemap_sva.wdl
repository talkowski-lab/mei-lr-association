version 1.0

import "utils/concat_files.wdl" as ConcatFiles
import "utils/download_avtable.wdl" as AVTableUtils
import "https://raw.githubusercontent.com/AoU-Multiomics-Analysis/susieR/refs/heads/main/workflows/susieRonly.wdl" as susieRonly


workflow FinemapSVAMultiomics {
  input {
    File SVAGTBed

    File eQTLAssocSummary
    File sQTLAssocSummary
    File pQTLAssocSummary

    File eQTL_SampleList
    File sQTL_SampleList
    File pQTL_SampleList

    File eQTL_Covariates
    File sQTL_Covariates
    File pQTL_Covariates

    File eQTL_TensorQTLPerm

    String eQTL_AVTable
    String sQTL_AVTable
    String pQTL_AVTable

    String WorkspaceNamespace
    String WorkspaceName

    String? Prefix

    Int CisDistance = 1000000
    Int Memory = 64
    Int NumPreempt = 3
    Float? MAF
  }

  String eQTL_FM_Prefix = if defined(Prefix) then select_first([Prefix]) + "_eQTL" else "eQTL"

  call SubsetSignificantSVAGT as eQTL_signifSVAs {
    input:
      SignatureType = "expression",
      QTLAssocSummary = eQTLAssocSummary,
      SVAGTBed = SVAGTBed,
      Prefix = eQTL_FM_Prefix
  }

  call AVTableUtils.DownloadAVTable as eQTL_DownloadAVTable {
    input:
      AVTableName = eQTL_AVTable,
      WorkspaceNamespace = WorkspaceNamespace,
      WorkspaceName = WorkspaceName
  }

  scatter (batch_idx in range(length(eQTL_signifSVAs.FeatureNames))) {
    String eQTL_FeatureName = eQTL_signifSVAs.FeatureNames[batch_idx]

    call GatherSNVFM_Input as eQTL_GetSNVFM_Input {
      input:
        AVTableFile = eQTL_DownloadAVTable.AVTableData,
        QTL_AVTable = eQTL_AVTable,
        GeneName = eQTL_FeatureName
    }

    call susieRonly.susieR as eQTL_Susie {
      input:
        GenotypeDosages = eQTL_GetSNVFM_Input.GenotypeDosages,
        GenotypeDosageIndex = eQTL_GetSNVFM_Input.GenotypeDosagesIndex,
        QTLCovariates = eQTL_Covariates,
        TensorQTLPermutations = eQTL_TensorQTLPerm,
        SampleList = eQTL_SampleList,
        PhenotypeBed = eQTL_GetSNVFM_Input.PhenotypeBed,
        CisDistance = CisDistance,
        AdditionalGenotypesBed = eQTL_signifSVAs.SignificantSVAGT,
        OutputPrefix = eQTL_FM_Prefix + "_" + eQTL_FeatureName,
        memory = Memory,
        NumPrempt = NumPreempt,
        MAF = MAF
    }
  }

  call ConcatFiles.ConcatParquet as eQTL_SusieParq_Concat {
    input:
      InputFiles = eQTL_Susie.SusieParquet,
      OutputName = eQTL_FM_Prefix + "_susie.parquet"
  }

  call ConcatFiles.ConcatParquet as eQTL_FullSusieParq_Concat {
    input:
      InputFiles = eQTL_Susie.FullSusieParquet,
      OutputName = eQTL_FM_Prefix + "_full_susie.parquet"
  }

  output {
    File eQTL_SusieParquet = eQTL_SusieParq_Concat.ConcatenatedFile
    File eQTL_FullSusieParquet = eQTL_FullSusieParq_Concat.ConcatenatedFile
    File eQTL_SignificantPairs = eQTL_signifSVAs.SignificantPairs
    File eQTL_SignificantSVAGT = eQTL_signifSVAs.SignificantSVAGT
  }
}


# Filter the Assoc Table to those with a pvalue below PValueThreshold, then
# join each significant (sva_id, feature) pair against the SVA genotype/
# repeat-length BED to produce the list of significant pairs plus one wide
# row per (variant, feature) combination.
# If SignatureType is splice, the leafcutter feature is collapsed to its gene id.
task SubsetSignificantSVAGT {
  input {
    String SignatureType
    File QTLAssocSummary
    File SVAGTBed
    Float PValueThreshold = 0.01
    String Prefix = "sig_sva"
    String ImageTag = "latest"
    Int MemoryGB = 4
    Int? DiskGB
  }

  Int auto_disk_size = ceil(size([QTLAssocSummary, SVAGTBed], "GB") * 2) + 10

  command <<<
    set -euo pipefail

    Rscript /scripts/subset_significant_sva_gt.R \
        --qtl-assoc-summary ~{QTLAssocSummary} \
        --sva-gt-bed ~{SVAGTBed} \
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

# Takes a string representing the table on terra and a string representing the gene name identifier and pulls the requested columns
task GatherSNVFM_Input {

  input {
    File AVTableFile
    String QTL_AVTable
    String GeneName
    String ImageTag = "latest"
    Int MemoryGB = 4
    Int? DiskGB
  }

  command <<<
    set -euo pipefail

    Rscript /scripts/gather_snv_fm_input.R \
        --avtable-file ~{AVTableFile} \
        --qtl-avtable ~{QTL_AVTable} \
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
