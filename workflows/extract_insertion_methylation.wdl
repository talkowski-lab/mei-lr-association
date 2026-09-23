version 1.0

workflow ExtractInsertionMethylation {
    input {
        File Bam
        File BamIndex
        File LociBed
        String Prefix
        String ImageTag = "latest"
        Float LenTolerancePerc = 0.05
        Int AnchorPad = 50
        String ModCodes = "m,h"
        Float ProbThreshold = 0.5
    }

    call ExtractMethylation {
        input:
            Bam = Bam,
            BamIndex = BamIndex,
            LociBed = LociBed,
            Prefix = Prefix,
            ImageTag = ImageTag,
            LenTolerancePerc = LenTolerancePerc,
            AnchorPad = AnchorPad,
            ModCodes = ModCodes,
            ProbThreshold = ProbThreshold
    }

    output {
        File MethylationTable = ExtractMethylation.MethylationTable
    }
}

task ExtractMethylation {
    input {
        File Bam
        File BamIndex
        File LociBed
        String Prefix
        String ImageTag = "latest"
        Float LenTolerancePerc = 0.05
        Int AnchorPad = 50
        String ModCodes = "m,h"
        Float ProbThreshold = 0.5
    }

    Int disk_size = ceil(size(Bam, "GB") * 1.5) + 20

    command <<<
        set -euo pipefail

        # Co-locate the index next to the BAM so pysam can find it.
        ln -s ~{Bam} input.bam
        ln -s ~{BamIndex} input.bam.bai

        python3 /scripts/extract_insertion_methylation.py \
            --bam input.bam \
            --loci-bed ~{LociBed} \
            --prefix ~{Prefix} \
            --len-tolerance-perc ~{LenTolerancePerc} \
            --anchor-pad ~{AnchorPad} \
            --mod-codes ~{ModCodes} \
            --prob-threshold ~{ProbThreshold}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-python_general:" + ImageTag
        memory: "16G"
        cpu: 4
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File MethylationTable = Prefix + ".insertion_methylation.tsv"
    }
}
