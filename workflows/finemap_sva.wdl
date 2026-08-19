version 1.0

import "finemap_sva_eqtl.wdl" as EQTL
import "finemap_sva_sqtl.wdl" as SQTL
import "finemap_sva_pqtl.wdl" as PQTL

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

  call EQTL.FinemapSVAeQTL as eQTL_Run {
    input:
      SVAGTBed = SVAGTBed,
      AssocSummary = eQTLAssocSummary,
      AVTable = eQTL_AVTable,
      SampleList = eQTL_SampleList,
      Covariates = eQTL_Covariates,
      TensorQTLPerm = eQTL_TensorQTLPerm,
      PhenotypeBed = eQTL_PhenotypeBed,
      Prefix = Prefix,
      PValueThreshold = PValueThreshold,
      SubsetImageTag = SubsetImageTag,
      CisDistance = CisDistance,
      Memory = Memory,
      NumPreempt = NumPreempt,
      MAF = MAF
  }

  call SQTL.FinemapSVAsQTL as sQTL_Run {
    input:
      SVAGTBed = SVAGTBed,
      AssocSummary = sQTLAssocSummary,
      AVTable = sQTL_AVTable,
      SampleList = sQTL_SampleList,
      Covariates = sQTL_Covariates,
      Prefix = Prefix,
      PValueThreshold = PValueThreshold,
      SubsetImageTag = SubsetImageTag,
      CisDistance = CisDistance,
      Memory = Memory,
      NumPreempt = NumPreempt,
      MAF = MAF
  }

  call PQTL.FinemapSVApQTL as pQTL_Run {
    input:
      SVAGTBed = SVAGTBed,
      AssocSummary = pQTLAssocSummary,
      AVTable = pQTL_AVTable,
      SampleList = pQTL_SampleList,
      Covariates = pQTL_Covariates,
      Prefix = Prefix,
      PValueThreshold = PValueThreshold,
      SubsetImageTag = SubsetImageTag,
      CisDistance = CisDistance,
      Memory = Memory,
      NumPreempt = NumPreempt,
      MAF = MAF
  }

  output {
    File eQTL_SusieParquet = eQTL_Run.SusieParquet
    Map[String, Array[File]] eQTL_SusieAuxFiles = eQTL_Run.SusieAuxFiles
    File sQTL_SusieParquet = sQTL_Run.SusieParquet
    Map[String, Array[File]] sQTL_SusieAuxFiles = sQTL_Run.SusieAuxFiles
    File pQTL_SusieParquet = pQTL_Run.SusieParquet
    Map[String, Array[File]] pQTL_SusieAuxFiles = pQTL_Run.SusieAuxFiles
  }
}
