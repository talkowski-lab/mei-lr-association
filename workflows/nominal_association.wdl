version 1.0

## Runs nominal cis-association between SVA repeat-length/genotype variants
## (a ProcessMEISVAData per-individual matrix, shared across all three QTL
## types) and each of the three molecular QTL feature sets -- eQTL
## (expression), sQTL (splicing), pQTL (protein) -- via
## scripts/nominal_association.R.

workflow NominalAssociation {
    input {
        File SVAFile
        File EQTLFeatureFile
        File? EQTLCovariateFile
        File SQTLFeatureFile
        File? SQTLCovariateFile
        File PQTLFeatureFile
        File? PQTLCovariateFile
        Int Window = 1000000
        String? Prefix
        String ImageTag = "latest"
    }

    String EQTLPrefix = if defined(Prefix) then select_first([Prefix]) + "_eQTL" else "eQTL"
    String SQTLPrefix = if defined(Prefix) then select_first([Prefix]) + "_sQTL" else "sQTL"
    String PQTLPrefix = if defined(Prefix) then select_first([Prefix]) + "_pQTL" else "pQTL"

    call NominalAssociationTask as EQTL {
        input:
            SVAFile = SVAFile,
            FeatureFile = EQTLFeatureFile,
            CovariateFile = EQTLCovariateFile,
            Window = Window,
            Prefix = EQTLPrefix,
            ImageTag = ImageTag
    }

    call NominalAssociationTask as SQTL {
        input:
            SVAFile = SVAFile,
            FeatureFile = SQTLFeatureFile,
            CovariateFile = SQTLCovariateFile,
            Window = Window,
            Prefix = SQTLPrefix,
            ImageTag = ImageTag
    }

    call NominalAssociationTask as PQTL {
        input:
            SVAFile = SVAFile,
            FeatureFile = PQTLFeatureFile,
            CovariateFile = PQTLCovariateFile,
            Window = Window,
            Prefix = PQTLPrefix,
            ImageTag = ImageTag
    }

    output {
        Array[File] NominalAssocModels = [EQTL.NominalAssocModels, SQTL.NominalAssocModels, PQTL.NominalAssocModels]
        File EQTLSummary = EQTL.NominalAssocSummary
        File SQTLSummary = SQTL.NominalAssocSummary
        File PQTLSummary = PQTL.NominalAssocSummary
    }
}

task NominalAssociationTask {
    input {
        File SVAFile
        File FeatureFile
        File? CovariateFile
        Int Window = 1000000
        String Prefix
        String ImageTag = "latest"
        Int MemoryGB = 8
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([SVAFile, FeatureFile, CovariateFile], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/nominal_association.R \
            --feature-file ~{FeatureFile} \
            --sva-file ~{SVAFile} \
            ~{"--covariate-file " + CovariateFile} \
            --window ~{Window} \
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
        File NominalAssocModels = "~{Prefix}_nominal_assoc.Rds"
        File NominalAssocSummary = "~{Prefix}_nominal_assoc_summary.txt.gz"
    }
}
