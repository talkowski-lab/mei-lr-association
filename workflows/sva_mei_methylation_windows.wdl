version 1.0

import "combine_and_window_variants.wdl" as CombineAndWindowVariants
import "bedtools_intersect.wdl" as BedtoolsIntersect
import "process_methylation_window.wdl" as ProcessMethylationWindow
import "utils/concat_files.wdl" as ConcatFiles

## End-to-end pipeline from a reference-SVA bed and an MEI bed, plus one
## individual's phased methylation beds (hap1/hap2), to a 100bp-binned
## methylation-by-distance-from-variant summary:
##
## 1. CombineAndWindowVariants: concatenates SVABed+MEIBed and flanks each
##    variant by FlankSize on each side (clamped to chromosome bounds),
##    carrying the original coordinates and strand through as extra columns.
## 2. BedtoolsIntersect (called once per haplotype): intersects each of
##    Hap1MethylBed/Hap2MethylBed against the windowed variant bed (`-wa -wb`
##    so both sides' columns are kept).
## 3. ProcessMethylationWindow (called once per haplotype, independently):
##    bins each haplotype's intersect bed into 100bp tiles by signed
##    distance from its variant and summarizes per (variant, haplotype,
##    distance bin) -- haplotype grouping comes from the methylation bed's
##    own embedded column, not from which file it was. Summarizing per-
##    haplotype before concatenating (rather than concatenating the raw
##    intersect beds first) keeps ProcessMethylationWindow decoupled from how
##    many haplotype files there are -- it can equally be run standalone on a
##    single combined/unphased methylation bed later.
## 4. ConcatFiles.ConcatenateDelim: concatenates the two per-haplotype
##    summary TSVs (each already has a header and a `haplotype` column) into
##    one final table.

workflow SVAMEIMethylationWindows {
    input {
        File SVABed
        File MEIBed
        File ChromSizes
        File Hap1MethylBed
        File Hap2MethylBed
        String Individual
        Int FlankSize = 10000
        Int BinSize = 100
        String ImageTag = "latest"
        String BedtoolsDockerImage = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
    }

    call CombineAndWindowVariants.CombineAndWindowVariants as CombineAndWindow {
        input:
            SVABed = SVABed,
            MEIBed = MEIBed,
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
