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
        ## Needed for the source-flank search only when Reference is a .mmi or absent
        File? ReferenceFasta
        File? ReferenceFastaIndex
        File? TransductionPaf
        File? RefRmsk
        Array[File] KnownSources = []
        File? Genes
        Boolean FlankSearch = true
        String Prefix = "mei"
        String ExtraArgs = ""
        String ImageTag = "latest"
        ## gs:// folder every output is also copied into; omit to skip
        String? OutputDir
    }

    call CatalogStructures {
        input:
            LimeaidTable = LimeaidTable,
            MeiTable = MeiTable,
            GenotypeTable = GenotypeTable,
            MeiVariantCol = MeiVariantCol,
            Reference = Reference,
            ReferenceFasta = ReferenceFasta,
            ReferenceFastaIndex = ReferenceFastaIndex,
            FlankSearch = FlankSearch,
            TransductionPaf = TransductionPaf,
            RefRmsk = RefRmsk,
            KnownSources = KnownSources,
            Genes = Genes,
            Prefix = Prefix,
            ExtraArgs = ExtraArgs,
            ImageTag = ImageTag
    }

    if (defined(OutputDir)) {
        call CopyOutputs {
            input:
                Files = select_all([
                    CatalogStructures.PerSequence,
                    CatalogStructures.PerVariant,
                    CatalogStructures.TransductionsFasta,
                    CatalogStructures.TransductionsPaf,
                    CatalogStructures.TransductionSources,
                    CatalogStructures.SourceSummary,
                    CatalogStructures.Transductions,
                    CatalogStructures.TransductionSegments
                ]),
                OutputDir = select_first([OutputDir]),
                ImageTag = ImageTag
        }
    }

    output {
        File PerSequence = CatalogStructures.PerSequence
        File PerVariant = CatalogStructures.PerVariant
        File TransductionsFasta = CatalogStructures.TransductionsFasta
        File? TransductionsPaf = CatalogStructures.TransductionsPaf
        File? TransductionSources = CatalogStructures.TransductionSources
        File? SourceSummary = CatalogStructures.SourceSummary
        File? Transductions = CatalogStructures.Transductions
        File? TransductionSegments = CatalogStructures.TransductionSegments
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
        ## FASTA for the source-flank search; defaults to Reference unless that is a .mmi
        File? ReferenceFasta
        File? ReferenceFastaIndex
        File? TransductionPaf
        File? RefRmsk
        Array[File] KnownSources = []
        File? Genes

        ## Second pass: BLAST untraced stretches >= MinSegment bp against the flanks of candidate source elements
        Boolean FlankSearch = true
        Int MinSegment = 20
        String FlankFamilies = "SVA,L1"
        String FlankMaxEvalue = "1e-5"
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
                            + size(select_all([Reference, ReferenceFasta, TransductionPaf, RefRmsk, Genes]), "GB")
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
            ~{"--reference-fasta " + ReferenceFasta} \
            ~{"--reference-fasta-index " + ReferenceFastaIndex} \
            --flank-search ~{FlankSearch} \
            --min-segment ~{MinSegment} \
            --flank-families '~{FlankFamilies}' \
            --flank-max-evalue ~{FlankMaxEvalue} \
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
        File? Transductions = "~{Prefix}.structure_transductions.tsv"
        File? TransductionSegments = "~{Prefix}.structure_transduction_segments.tsv"
    }
}

task CopyOutputs {
    input {
        Array[File] Files
        String OutputDir
        String ImageTag = "latest"
    }

    command <<<
        set -euo pipefail
        gcloud storage cp ~{sep=" " Files} "~{sub(OutputDir, "/+$", "")}/"
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
        memory: "2 GB"
        cpu: 1
        disks: "local-disk " + (ceil(size(Files, "GB")) + 10) + " HDD"
        preemptible: 2
        maxRetries: 1
    }
}
