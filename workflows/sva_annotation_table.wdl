version 1.0

## Builds the per-SVA annotation table (scripts/sva_annotation_table.R) that
## notebooks/sva_enrichment_permutation.Rmd reads: one row per polymorphic and
## fixed SVA with allele count/frequency, repeat-length statistics, structural
## class from MEIStructuralCatalog, and every genomic annotation's value, plus
## the annotation bundle (.rds) the notebook's positional null draws from.
##
## Every input is optional. One left unset falls back to the script's default:
## the file the notebook used, read with `gcloud storage cp` from GcloudBucket
## (bucket-relative paths) or RefFolder within it (reference files), with public
## UCSC / GENCODE copies for the genome, FASTA and GTF.

workflow SVAAnnotationTable {
    input {
        ## Workspace bucket the default inputs are read from (with or without gs://)
        String? GcloudBucket
        String RefFolder = "ref"
        ## Clustered polymorphic SVA BED (chrom, start, end, REF, ALT, ID, MEI_type, MEI_subfamily, length, orientation, intact_status, cluster_id)
        File? MeiBed
        ## ProcessMEISVAData SVA_MEI_lengths_tidy.txt.gz
        File? MeiLengths
        ## SVA_regions_van_bree.bed (or _flanked.txt)
        File? RefTable
        ## ProcessRefSVAData ref_SVA_length_filtered.txt.gz
        File? RefLengths
        ## FilterMEIVCF SVA_gt_tidy_clustered.txt.gz
        File? GenotypeTable
        ## participant, sex: chrX/chrY allele numbers and sex filtering
        File? SampleMetadata
        ## MEIStructuralCatalog <prefix>.structure_per_variant.tsv
        File? CatalogPerVariant
        ## Polymorphic Alu / L1, for their own tests: each family's clustered bed (same
        ## 12 columns as MeiBed) and FilterMEIVCF's <MEIType>_gt_tidy_clustered.txt.gz
        File? AluBed
        File? AluGenotypeTable
        File? L1Bed
        File? L1GenotypeTable
        File? Genome
        ## GC from bedtools nuc (uncompressed FASTA); or give GcBedgraph
        File? Reference
        File? GcBedgraph
        File? GencodeGtf
        ## gnomAD v4.1 constraint metrics, for the gene_constraint annotation
        File? LoeufTable
        File? ExcludeBed
        File? IncludeBed
        ## hg38 OK-seq RFD; without it the published HeLa hg19 track is lifted over
        File? RfdBedgraph
        File? OkseqWatson
        File? OkseqCrick
        ## Prebuilt replication timing; without it RT is built from 4DN Repli-seq
        File? RepliseqBedgraph
        String Prefix = "sva"
        ## Any other scripts/sva_annotation_table.R option, e.g. "--flanks 10000,2000 --rfd-flip true"
        String ExtraArgs = ""
        String ImageTag = "latest"
        Int MemoryGB = 32
        Int DiskGB = 100
        ## gs:// folder every output is also copied into; omit to skip
        String? OutputDir
    }

    call BuildAnnotationTable {
        input:
            GcloudBucket = GcloudBucket,
            RefFolder = RefFolder,
            MeiBed = MeiBed,
            MeiLengths = MeiLengths,
            RefTable = RefTable,
            RefLengths = RefLengths,
            GenotypeTable = GenotypeTable,
            SampleMetadata = SampleMetadata,
            CatalogPerVariant = CatalogPerVariant,
            AluBed = AluBed,
            AluGenotypeTable = AluGenotypeTable,
            L1Bed = L1Bed,
            L1GenotypeTable = L1GenotypeTable,
            Genome = Genome,
            Reference = Reference,
            GcBedgraph = GcBedgraph,
            GencodeGtf = GencodeGtf,
            LoeufTable = LoeufTable,
            ExcludeBed = ExcludeBed,
            IncludeBed = IncludeBed,
            RfdBedgraph = RfdBedgraph,
            OkseqWatson = OkseqWatson,
            OkseqCrick = OkseqCrick,
            RepliseqBedgraph = RepliseqBedgraph,
            Prefix = Prefix,
            ExtraArgs = ExtraArgs,
            ImageTag = ImageTag,
            MemoryGB = MemoryGB,
            DiskGB = DiskGB
    }

    if (defined(OutputDir)) {
        call CopyOutputs {
            input:
                Files = flatten([[BuildAnnotationTable.AnnotationTable, BuildAnnotationTable.AnnotationBundle],
                                 select_all([BuildAnnotationTable.AluAnnotationTable, BuildAnnotationTable.L1AnnotationTable]),
                                 BuildAnnotationTable.Tracks]),
                OutputDir = select_first([OutputDir]),
                ImageTag = ImageTag
        }
    }

    output {
        File AnnotationTable = BuildAnnotationTable.AnnotationTable
        File AnnotationBundle = BuildAnnotationTable.AnnotationBundle
        File? AluAnnotationTable = BuildAnnotationTable.AluAnnotationTable
        File? L1AnnotationTable = BuildAnnotationTable.L1AnnotationTable
        Array[File] Tracks = BuildAnnotationTable.Tracks
    }
}

task BuildAnnotationTable {
    input {
        String? GcloudBucket
        String RefFolder = "ref"
        File? MeiBed
        File? MeiLengths
        File? RefTable
        File? RefLengths
        File? GenotypeTable
        File? SampleMetadata
        File? CatalogPerVariant
        File? AluBed
        File? AluGenotypeTable
        File? L1Bed
        File? L1GenotypeTable
        File? Genome
        File? Reference
        File? GcBedgraph
        File? GencodeGtf
        File? LoeufTable
        File? ExcludeBed
        File? IncludeBed
        File? RfdBedgraph
        File? OkseqWatson
        File? OkseqCrick
        File? RepliseqBedgraph
        String Prefix = "sva"
        String ExtraArgs = ""
        String ImageTag = "latest"
        Int MemoryGB = 32
        ## default inputs are downloaded at run time (hg38.fa alone is ~3 GB, decompressed if gzipped)
        Int DiskGB = 100
    }

    command <<<
        set -euo pipefail

        Rscript /scripts/sva_annotation_table.R \
            ~{"--gcloud-bucket " + GcloudBucket} \
            --ref-folder '~{RefFolder}' \
            ~{"--mei-bed " + MeiBed} \
            ~{"--mei-lengths " + MeiLengths} \
            ~{"--ref-table " + RefTable} \
            ~{"--ref-lengths " + RefLengths} \
            ~{"--genotype-table " + GenotypeTable} \
            ~{"--genome " + Genome} \
            ~{"--sample-metadata " + SampleMetadata} \
            ~{"--catalog-per-variant " + CatalogPerVariant} \
            ~{"--alu-bed " + AluBed} \
            ~{"--alu-genotype-table " + AluGenotypeTable} \
            ~{"--l1-bed " + L1Bed} \
            ~{"--l1-genotype-table " + L1GenotypeTable} \
            ~{"--reference " + Reference} \
            ~{"--gc-bedgraph " + GcBedgraph} \
            ~{"--gencode-gtf " + GencodeGtf} \
            ~{"--loeuf-table " + LoeufTable} \
            ~{"--exclude-bed " + ExcludeBed} \
            ~{"--include-bed " + IncludeBed} \
            ~{"--rfd-bedgraph " + RfdBedgraph} \
            ~{"--okseq-watson " + OkseqWatson} \
            ~{"--okseq-crick " + OkseqCrick} \
            ~{"--repliseq-bigwig " + RepliseqBedgraph} \
            --outdir . \
            --prefix ~{Prefix} \
            ~{ExtraArgs}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: 2
        disks: "local-disk " + DiskGB + " SSD"
        preemptible: 1
        maxRetries: 1
    }

    output {
        File AnnotationTable = "~{Prefix}.sva_annotation_table.tsv.gz"
        File AnnotationBundle = "~{Prefix}.sva_annotations.rds"
        File? AluAnnotationTable = "~{Prefix}.alu_annotation_table.tsv.gz"
        File? L1AnnotationTable = "~{Prefix}.l1_annotation_table.tsv.gz"
        Array[File] Tracks = glob("~{Prefix}.*.bed*")
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
