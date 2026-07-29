version 1.0

## Aligns both haplotype assemblies of a diploid assembly to a reference with
## minimap2, producing one coordinate-sorted, indexed BAM per haplotype. The
## two haplotypes are aligned independently (and run concurrently under
## Cromwell) since they don't depend on each other.

workflow AlignDiploidAssemblyToRef {
    input {
        File Hap1Contigs
        File Hap2Contigs
        File Reference
        String Prefix
        String MinimapPreset = "asm5"
        String ImageTag = "latest"
    }

    call AlignAssembly as AlignHap1 {
        input:
            Contigs = Hap1Contigs,
            Reference = Reference,
            Prefix = Prefix + ".hap1",
            MinimapPreset = MinimapPreset,
            ImageTag = ImageTag
    }

    call AlignAssembly as AlignHap2 {
        input:
            Contigs = Hap2Contigs,
            Reference = Reference,
            Prefix = Prefix + ".hap2",
            MinimapPreset = MinimapPreset,
            ImageTag = ImageTag
    }

    output {
        File Hap1Bam = AlignHap1.Bam
        File Hap1BamIndex = AlignHap1.BamIndex
        File Hap2Bam = AlignHap2.Bam
        File Hap2BamIndex = AlignHap2.BamIndex
    }
}

task AlignAssembly {
    input {
        File Contigs
        File Reference
        String Prefix
        String MinimapPreset = "asm5"
        String ImageTag = "latest"
        Int CPU = 8
        Int MemoryGB = 16
        Int? DiskGB
    }

    # +2x Contigs size for the intermediate SAM file (uncompressed, written to disk rather
    # than piped straight into samtools sort -- see comment in command below).
    Int auto_disk_size = ceil((size(Contigs, "GB") * 3 + size(Reference, "GB")) * 4) + 20

    command <<<
        set -euo pipefail

        # minimap2's SAM output is written to a file rather than piped directly into
        # samtools sort: piping the two together can race on writing/reading the SAM
        # header when the reference has many contigs (e.g. hg38-style chrUn_*/alt/random
        # scaffolds), corrupting one @SQ line and making samtools sort fail with
        # "bad or missing LN tag". Writing to a file first fully decouples the two steps.
        minimap2 -ax ~{MinimapPreset} -t ~{CPU} ~{Reference} ~{Contigs} > aligned.sam

        samtools sort -@ ~{CPU} -o ~{Prefix}.bam aligned.sam

        samtools index ~{Prefix}.bam
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-bioinformatics:" + ImageTag
        memory: MemoryGB + " GB"
        cpu: CPU
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File Bam = "~{Prefix}.bam"
        File BamIndex = "~{Prefix}.bam.bai"
    }
}
