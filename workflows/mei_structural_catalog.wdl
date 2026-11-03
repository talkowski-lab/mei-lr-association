version 1.0

## Catalogs the structural forms of non-reference Alu/L1/SVA insertions from
## L1ME-AID/INTACT_MEI tables via scripts/mei_structural_catalog.R, and, when
## Reference or TransductionPaf is given, traces the genomic sources of their
## transduced sequence. Analysis of the outputs:
## notebooks/mei_structural_catalog_analysis.Rmd.

workflow MEIStructuralCatalog {
    input {
        ## L1ME-AID/INTACT_MEI TSV covering every family
        File LimeaidTable
        ## Clustered MEI allele table: variant grouping and family (MEI_type)
        File MeiTable
        ## CHROM, POS, ALT_md5, indiv, GT, hap1, hap2, ...
        File GenotypeTable
        String MeiVariantCol = "var_id"
        File? Reference
        File? TransductionPaf
        File? RefRmsk
        Array[File] KnownSources = []
        File? Genes
        String Prefix = "mei"
        String ExtraArgs = ""
        String ImageTag = "latest"
    }

    call CatalogStructures {
        input:
            LimeaidTable = LimeaidTable,
            MeiTable = MeiTable,
            GenotypeTable = GenotypeTable,
            MeiVariantCol = MeiVariantCol,
            Reference = Reference,
            TransductionPaf = TransductionPaf,
            RefRmsk = RefRmsk,
            KnownSources = KnownSources,
            Genes = Genes,
            Prefix = Prefix,
            ExtraArgs = ExtraArgs,
            ImageTag = ImageTag
    }

    output {
        File PerSequence = CatalogStructures.PerSequence
        File PerVariant = CatalogStructures.PerVariant
        File TransductionsFasta = CatalogStructures.TransductionsFasta
        File? TransductionsPaf = CatalogStructures.TransductionsPaf
        File? TransductionSources = CatalogStructures.TransductionSources
        File? SourceSummary = CatalogStructures.SourceSummary
    }
}

task CatalogStructures {
    input {
        File LimeaidTable
        File MeiTable
        File GenotypeTable
        String MeiVariantCol = "var_id"
        ## AF denominator is 2 x this; 0 = distinct indiv in GenotypeTable
        Int NIndividuals = 0
        ## FASTA or minimap2 .mmi; a FASTA is indexed on the fly, which needs ~16 GB for GRCh38
        File? Reference
        File? TransductionPaf
        File? RefRmsk
        Array[File] KnownSources = []
        File? Genes

        String Max5pMissing = "Alu=20,L1=50,SVA=50"
        String Max3pMissing = "Alu=20,L1=50,SVA=50"
        Int MinTransduction = 30
        Int MinPolyA = 8
        Int TerminalSlack = 30
        Int MinInversion = 20
        String MinimapPreset = "sr"
        ## Any other scripts/mei_structural_catalog.R option, e.g. "--source-window 5000 --run-ava false"
        String ExtraArgs = ""

        String Prefix = "mei"
        String ImageTag = "latest"
        Int CPU = 4
        Int MemoryGB = 32
        Int? DiskGB
    }

    Int auto_disk_gb = ceil(size(LimeaidTable, "GB") + size(MeiTable, "GB") + size(GenotypeTable, "GB")
                            + size(select_all([Reference, TransductionPaf, RefRmsk, Genes]), "GB")
                            + size(KnownSources, "GB")) * 2 + 20

    command <<<
        set -euo pipefail

        KNOWN_SOURCES='~{sep="," KnownSources}'

        Rscript /scripts/mei_structural_catalog.R \
            --limeaid-table ~{LimeaidTable} \
            --mei-table ~{MeiTable} \
            --genotype-table ~{GenotypeTable} \
            --mei-variant-col ~{MeiVariantCol} \
            --n-individuals ~{NIndividuals} \
            ~{"--reference " + Reference} \
            ~{"--transduction-paf " + TransductionPaf} \
            ~{"--ref-rmsk " + RefRmsk} \
            ${KNOWN_SOURCES:+--known-sources "$KNOWN_SOURCES"} \
            ~{"--genes " + Genes} \
            --max-5p-missing '~{Max5pMissing}' \
            --max-3p-missing '~{Max3pMissing}' \
            --min-transduction ~{MinTransduction} \
            --min-polya ~{MinPolyA} \
            --terminal-slack ~{TerminalSlack} \
            --min-inversion ~{MinInversion} \
            --minimap2-preset ~{MinimapPreset} \
            --threads ~{CPU} \
            --outdir . \
            --prefix ~{Prefix} \
            ~{ExtraArgs}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: CPU
        disks: "local-disk " + select_first([DiskGB, auto_disk_gb]) + " SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File PerSequence = "~{Prefix}.structure_per_sequence.tsv.gz"
        File PerVariant = "~{Prefix}.structure_per_variant.tsv"
        File TransductionsFasta = "~{Prefix}.structure_transductions.fa"
        File? TransductionsPaf = "~{Prefix}.structure_transductions.paf"
        File? TransductionSources = "~{Prefix}.structure_transduction_sources.tsv"
        File? SourceSummary = "~{Prefix}.structure_source_summary.tsv"
    }
}
