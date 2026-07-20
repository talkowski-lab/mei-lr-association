version 1.0

workflow ExtractInsertionMethylation {
    input {
        File Bam
        File BamIndex
        File LociBed
        String Branch = "main"
        Int MinInsLen = 2000
        Int MaxInsLen = 4000
        Int AnchorPad = 100
        String ModCodes = "m,h"
        Float ProbThreshold = 0.5
    }

    call ExtractMethylation {
        input:
            Bam = Bam,
            BamIndex = BamIndex,
            LociBed = LociBed,
            Branch = Branch,
            MinInsLen = MinInsLen,
            MaxInsLen = MaxInsLen,
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
        String Branch = "main"
        Int MinInsLen = 2000
        Int MaxInsLen = 4000
        Int AnchorPad = 100
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
            --output insertion_methylation.tsv \
            --min-ins-len ~{MinInsLen} \
            --max-ins-len ~{MaxInsLen} \
            --anchor-pad ~{AnchorPad} \
            --mod-codes ~{ModCodes} \
            --prob-threshold ~{ProbThreshold}
    >>>

    runtime {
        docker: "ghcr.io/alyenkin/insertion-methylation:" + Branch
        memory: "16G"
        cpu: 4
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File MethylationTable = "insertion_methylation.tsv"
    }
}
