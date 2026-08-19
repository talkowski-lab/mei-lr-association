version 1.0

import "concat_files_with_script.wdl" as ConcatFiles
import "utils/build_file_map.wdl" as FileMapUtils
import "utils/finemap_sva_tasks.wdl" as FinemapTasks

workflow FinemapSVAeQTL {
  input {
    File SVAGTBed
    File AssocSummary
    File AVTable
    File SampleList
    File Covariates
    File TensorQTLPerm
    File PhenotypeBed

    String? Prefix
    Float? PValueThreshold
    String? SubsetImageTag
    Int CisDistance = 1000000
    Int Memory = 64
    Int NumPreempt = 3
    Float? MAF
  }

  String FM_Prefix = if defined(Prefix) then select_first([Prefix]) + "_eQTL" else "eQTL"

  call FinemapTasks.SubsetSignificantSVAGT as signifSVAs {
    input:
      SignatureType = "expression",
      QTLAssocSummary = AssocSummary,
      SVAGTBed = SVAGTBed,
      AVTable = AVTable,
      PValueThreshold = PValueThreshold,
      Prefix = FM_Prefix,
      ImageTag = SubsetImageTag
  }

  scatter (batch_idx in range(length(signifSVAs.FeatureNames))) {
    String FeatureName = signifSVAs.FeatureNames[batch_idx]

    call FinemapTasks.GatherSNVFM_Input as GetSNVFM_Input {
      input:
        AVTableFile = AVTable,
        GeneName = FeatureName
    }

    call FinemapTasks.susieR as Susie {
      input:
        GenotypeDosages = GetSNVFM_Input.GenotypeDosages,
        GenotypeDosageIndex = GetSNVFM_Input.GenotypeDosagesIndex,
        QTLCovariates = Covariates,
        TensorQTLPermutations = TensorQTLPerm,
        SampleList = SampleList,
        PhenotypeBed = PhenotypeBed,
        CisDistance = CisDistance,
        AdditionalGenotypesBed = signifSVAs.SignificantSVAGT,
        OutputPrefix = FeatureName,
        memory = Memory,
        NumPrempt = NumPreempt,
        MAF = MAF
    }
  }

  call ConcatFiles.ConcatenateAndProcessFiles as SusieParq_Concat {
    input:
      InputFiles = Susie.SusieParquet,
      OutputName = FM_Prefix + "_susie.parquet",
      FileType = "parquet",
      BatchSize = 100,
      BatchMemoryGB = 8,
      MergeMemoryGB = 32
  }

  call FileMapUtils.BuildFileMap as SusieAuxFileMap {
    input:
      Keys = FeatureName,
      Values = [Susie.lbfParquet, Susie.FullSusieParquet, Susie.VariantPositionSummary, Susie.SusieObject]
  }

  output {
    File SusieParquet = SusieParq_Concat.ConcatenatedFile
    Map[String, Array[File]] SusieAuxFiles = SusieAuxFileMap.FileMap
  }
}
