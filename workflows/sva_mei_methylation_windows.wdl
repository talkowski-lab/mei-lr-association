version 1.0

import "combine_and_window_variants.wdl" as CombineAndWindowVariants
import "bedtools_intersect.wdl" as BedtoolsIntersect
import "process_methylation_window.wdl" as ProcessMethylationWindow
import "utils/concat_files.wdl" as ConcatFiles

## End-to-end pipeline from a reference-SVA bed, an MEI bed, and a third bed
## (same format as the MEI bed), plus one individual's phased methylation
## beds (hap1/hap2), to a fixed-width-binned methylation-by-distance-from-
## variant summary: CombineAndWindowVariants -> BedtoolsIntersect (once per
## haplotype) -> ProcessMethylationWindow (once per haplotype) ->
## ConcatFiles.ConcatenateDelim.

workflow SVAMEIMethylationWindows {
    input {
        File SVABed
        File MEIBed
        File OtherBed
        File ChromSizes
        File Hap1MethylBed
        File Hap2MethylBed
        String Individual
        Int FlankSize = 5000
        Int BinSize = 200
        String ImageTag = "latest"
        String BedtoolsDockerImage = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
    }

    call CombineAndWindowVariants.CombineAndWindowVariants as CombineAndWindow {
        input:
            SVABed = SVABed,
            MEIBed = MEIBed,
            OtherBed = OtherBed,
            ChromSizes = ChromSizes,
            FlankSize = FlankSize,
            Prefix = Individual,
            DockerImage = BedtoolsDockerImage
    }

    call BedtoolsIntersect.BedtoolsIntersect as IntersectHap1 {
        input:
            FileA = Hap1MethylBed,
            FileB = CombineAndWindow.WindowedVariantBed,
            Params = "-wa -wb",
            Prefix = Individual + ".hap1",
            DockerImage = BedtoolsDockerImage
    }

    call BedtoolsIntersect.BedtoolsIntersect as IntersectHap2 {
        input:
            FileA = Hap2MethylBed,
            FileB = CombineAndWindow.WindowedVariantBed,
            Params = "-wa -wb",
            Prefix = Individual + ".hap2",
            DockerImage = BedtoolsDockerImage
    }

    call ProcessMethylationWindow.ProcessMethylationWindow as SummarizeHap1 {
        input:
            Individual = Individual,
            IntersectBed = IntersectHap1.IntersectOutput,
            BinSize = BinSize,
            FlankSize = FlankSize,
            Prefix = Individual + ".hap1",
            ImageTag = ImageTag
    }

    call ProcessMethylationWindow.ProcessMethylationWindow as SummarizeHap2 {
        input:
            Individual = Individual,
            IntersectBed = IntersectHap2.IntersectOutput,
            BinSize = BinSize,
            FlankSize = FlankSize,
            Prefix = Individual + ".hap2",
            ImageTag = ImageTag
    }

    call ConcatFiles.ConcatenateDelim as MergeHaplotypeSummaries {
        input:
            InputFiles = [SummarizeHap1.MethylationWindowSummary, SummarizeHap2.MethylationWindowSummary],
            HasHeader = true,
            OutputName = Individual + ".methylation_window_summary.tsv",
            ImageTag = ImageTag
    }

    output {
        File WindowedVariantBed = CombineAndWindow.WindowedVariantBed
        File Hap1MethylationWindowSummary = SummarizeHap1.MethylationWindowSummary
        File Hap2MethylationWindowSummary = SummarizeHap2.MethylationWindowSummary
        File MethylationWindowSummary = MergeHaplotypeSummaries.ConcatenatedFile
    }
}
