version 1.0

import "concat_files_with_script.wdl" as ConcatFiles
import "utils/build_file_map.wdl" as FileMapUtils
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
    File eQTL_PhenotypeBed

    # Already-downloaded QTL data tables (id column first).
    File eQTL_AVTable
    File sQTL_AVTable
    File pQTL_AVTable

    String? Prefix

    Float? PValueThreshold
    String? SubsetImageTag
    Int CisDistance = 1000000
    Int Memory = 64
    Int NumPreempt = 3
    Float? MAF
  }

  ## ----------------------------------------------------------------------------------------
  # eQTL Section
  String eQTL_FM_Prefix = if defined(Prefix) then select_first([Prefix]) + "_eQTL" else "eQTL"

  call SubsetSignificantSVAGT as eQTL_signifSVAs {
    input:
      SignatureType = "expression",
      QTLAssocSummary = eQTLAssocSummary,
      SVAGTBed = SVAGTBed,
      AVTable = eQTL_AVTable,
      PValueThreshold = PValueThreshold,
      Prefix = eQTL_FM_Prefix,
      ImageTag = SubsetImageTag
  }

  scatter (batch_idx in range(length(eQTL_signifSVAs.FeatureNames))) {
    String eQTL_FeatureName = eQTL_signifSVAs.FeatureNames[batch_idx]

    call GatherSNVFM_Input as eQTL_GetSNVFM_Input {
      input:
        AVTableFile = eQTL_AVTable,
        GeneName = eQTL_FeatureName
    }

    call susieRonly.susieR as eQTL_Susie {
      input:
        GenotypeDosages = eQTL_GetSNVFM_Input.GenotypeDosages,
        GenotypeDosageIndex = eQTL_GetSNVFM_Input.GenotypeDosagesIndex,
        QTLCovariates = eQTL_Covariates,
        TensorQTLPermutations = eQTL_TensorQTLPerm,
        SampleList = eQTL_SampleList,
        PhenotypeBed = eQTL_PhenotypeBed,
        CisDistance = CisDistance,
        AdditionalGenotypesBed = eQTL_signifSVAs.SignificantSVAGT,
        OutputPrefix = eQTL_FeatureName,
        memory = Memory,
        NumPrempt = NumPreempt,
        MAF = MAF
    }
  }

  call ConcatFiles.ConcatenateAndProcessFiles as eQTL_SusieParq_Concat {
    input:
      InputFiles = eQTL_Susie.SusieParquet,
      OutputName = eQTL_FM_Prefix + "_susie.parquet",
      FileType = "parquet",
      BatchSize = 100,
      BatchMemoryGB = 8,
      MergeMemoryGB = 32
  }

  call FileMapUtils.BuildFileMap as eQTL_SusieAuxFileMap {
    input:
      Keys = eQTL_FeatureName,
      Values = [eQTL_Susie.lbfParquet, eQTL_Susie.FullSusieParquet, eQTL_Susie.VariantPositionSummary, eQTL_Susie.SusieObject]
  }

  # call ConcatFiles.ConcatenateAndProcessFiles as eQTL_FullSusieParq_Concat {
  #   input:
  #     InputFiles = eQTL_Susie.FullSusieParquet,
  #     OutputName = eQTL_FM_Prefix + "_full_susie.parquet",
  #     FileType = "parquet",
  #     BatchSize = 100,
  #     BatchMemoryGB = 8,
  #     MergeMemoryGB = 32
  # }

  ## ----------------------------------------------------------------------------------------
  # sQTL Section
  String sQTL_FM_Prefix = if defined(Prefix) then select_first([Prefix]) + "_sQTL" else "sQTL"

  call SubsetSignificantSVAGT as sQTL_signifSVAs {
    input:
      SignatureType = "splice",
      QTLAssocSummary = sQTLAssocSummary,
      SVAGTBed = SVAGTBed,
      AVTable = sQTL_AVTable,
      PValueThreshold = PValueThreshold,
      Prefix = sQTL_FM_Prefix,
      ImageTag = SubsetImageTag
  }

  scatter (batch_idx in range(length(sQTL_signifSVAs.FeatureNames))) {
    String sQTL_FeatureName = sQTL_signifSVAs.FeatureNames[batch_idx]

    call GatherSNVFM_Input as sQTL_GetSNVFM_Input {
      input:
        AVTableFile = sQTL_AVTable,
        GeneName = sQTL_FeatureName
    }

    call susieRonly.susieR as sQTL_Susie {
      input:
        GenotypeDosages = sQTL_GetSNVFM_Input.GenotypeDosages,
        GenotypeDosageIndex = sQTL_GetSNVFM_Input.GenotypeDosagesIndex,
        QTLCovariates = sQTL_Covariates,
        TensorQTLPermutations = sQTL_GetSNVFM_Input.TensorQTLPermutations,
        SampleList = sQTL_SampleList,
        PhenotypeBed = sQTL_GetSNVFM_Input.PhenotypeBed,
        CisDistance = CisDistance,
        AdditionalGenotypesBed = sQTL_signifSVAs.SignificantSVAGT,
        OutputPrefix = sQTL_FeatureName,
        memory = Memory,
        NumPrempt = NumPreempt,
        MAF = MAF
    }
  }

  call ConcatFiles.ConcatenateAndProcessFiles as sQTL_SusieParq_Concat {
    input:
      InputFiles = sQTL_Susie.SusieParquet,
      OutputName = sQTL_FM_Prefix + "_susie.parquet",
      FileType = "parquet",
      BatchSize = 100,
      BatchMemoryGB = 8,
      MergeMemoryGB = 32
  }

  call FileMapUtils.BuildFileMap as sQTL_SusieAuxFileMap {
    input:
      Keys = sQTL_FeatureName,
      Values = [sQTL_Susie.lbfParquet, sQTL_Susie.FullSusieParquet, sQTL_Susie.VariantPositionSummary, sQTL_Susie.SusieObject]
  }

  ## ----------------------------------------------------------------------------------------
  # pQTL Section
  String pQTL_FM_Prefix = if defined(Prefix) then select_first([Prefix]) + "_pQTL" else "pQTL"

  call SubsetSignificantSVAGT as pQTL_signifSVAs {
    input:
      SignatureType = "protein",
      QTLAssocSummary = pQTLAssocSummary,
      SVAGTBed = SVAGTBed,
      AVTable = pQTL_AVTable,
      PValueThreshold = PValueThreshold,
      Prefix = pQTL_FM_Prefix,
      ImageTag = SubsetImageTag
  }

  scatter (batch_idx in range(length(pQTL_signifSVAs.FeatureNames))) {
    String pQTL_FeatureName = pQTL_signifSVAs.FeatureNames[batch_idx]

    call GatherSNVFM_Input as pQTL_GetSNVFM_Input {
      input:
        AVTableFile = pQTL_AVTable,
        GeneName = pQTL_FeatureName
    }

    call susieRonly.susieR as pQTL_Susie {
      input:
        GenotypeDosages = pQTL_GetSNVFM_Input.GenotypeDosages,
        GenotypeDosageIndex = pQTL_GetSNVFM_Input.GenotypeDosagesIndex,
        QTLCovariates = pQTL_Covariates,
        TensorQTLPermutations = pQTL_GetSNVFM_Input.TensorQTLPermutations,
        SampleList = pQTL_SampleList,
        PhenotypeBed = pQTL_GetSNVFM_Input.PhenotypeBed,
        CisDistance = CisDistance,
        AdditionalGenotypesBed = pQTL_signifSVAs.SignificantSVAGT,
        OutputPrefix = pQTL_FeatureName,
        memory = Memory,
        NumPrempt = NumPreempt,
        MAF = MAF
    }
  }

  call ConcatFiles.ConcatenateAndProcessFiles as pQTL_SusieParq_Concat {
    input:
      InputFiles = pQTL_Susie.SusieParquet,
      OutputName = pQTL_FM_Prefix + "_susie.parquet",
      FileType = "parquet",
      BatchSize = 100,
      BatchMemoryGB = 8,
      MergeMemoryGB = 32
  }

  call FileMapUtils.BuildFileMap as pQTL_SusieAuxFileMap {
    input:
      Keys = pQTL_FeatureName,
      Values = [pQTL_Susie.lbfParquet, pQTL_Susie.FullSusieParquet, pQTL_Susie.VariantPositionSummary, pQTL_Susie.SusieObject]
  }

  output {
    File eQTL_SusieParquet = eQTL_SusieParq_Concat.ConcatenatedFile
    Map[String, Array[File]] eQTL_SusieAuxFiles = eQTL_SusieAuxFileMap.FileMap
    File sQTL_SusieParquet = sQTL_SusieParq_Concat.ConcatenatedFile
    Map[String, Array[File]] sQTL_SusieAuxFiles = sQTL_SusieAuxFileMap.FileMap
    File pQTL_SusieParquet = pQTL_SusieParq_Concat.ConcatenatedFile
    Map[String, Array[File]] pQTL_SusieAuxFiles = pQTL_SusieAuxFileMap.FileMap
  }
}


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
