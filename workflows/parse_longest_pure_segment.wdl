version 1.0

## For each SVA (ID) x repeat region (hexamer/VNTR_1/VNTR_2/VNTR_3), finds the
## longest pure segment -- the longest run of consecutive repeat-unit calls
## that all match the same motif -- and reports that motif and its length,
## via scripts/parse_longest_pure_segment.R.

workflow ParseLongestPureSegment {
    input {
        File RepeatSeqTable
        String Prefix
        String ImageTag = "latest"
    }

    call FindLongestPureSegment {
        input:
            RepeatSeqTable = RepeatSeqTable,
            Prefix = Prefix,
            ImageTag = ImageTag
    }

    output {
        File LongestPureSegments = FindLongestPureSegment.LongestPureSegments
    }
}

task FindLongestPureSegment {
    input {
        File RepeatSeqTable
        String Prefix
        String ImageTag = "latest"
        Int MemoryGB = 4
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size(RepeatSeqTable, "GB") * 2) + 10

    command <<<
        set -euo pipefail

        Rscript /scripts/parse_longest_pure_segment.R \
            --repeat-seq-table ~{RepeatSeqTable} \
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
        File LongestPureSegments = "~{Prefix}_lps_length.txt"
    }
}
