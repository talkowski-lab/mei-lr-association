version 1.0

## Thin wrapper around `bedtools intersect`. Deliberately does not expose a
## dedicated input for every bedtools intersect flag -- Params is passed
## through to the command line as-is, so any combination of options (-v, -wa,
## -wb, -f, -r, -s, -sorted -g, etc.) is available without the WDL needing to
## track bedtools' option surface.

workflow BedtoolsIntersect {
    input {
        File FileA
        File? FileAIndex
        File FileB
        File? FileBIndex
        String Params = ""
        String Prefix = "intersect"
        String DockerImage = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1"
    }

    call Intersect {
        input:
            FileA = FileA,
            FileAIndex = FileAIndex,
            FileB = FileB,
            FileBIndex = FileBIndex,
            Params = Params,
            Prefix = Prefix,
            DockerImage = DockerImage
    }

    output {
        File IntersectOutput = Intersect.Output
    }
}

task Intersect {
    input {
        File FileA
        File? FileAIndex
        File FileB
        File? FileBIndex
        String Params = ""
        String Prefix = "intersect"
        String DockerImage
        Int? DiskGB
    }

    Int auto_disk_size = ceil(size([FileA, FileB], "GB") * 2) + 10

    command <<<
        set -euo pipefail

        # bedtools intersect streams -a/-b directly (bgzip-aware) and never seeks via
        # an index, so FileAIndex/FileBIndex aren't read here -- they're accepted and
        # symlinked purely so a bgzip+tabix pair can be passed straight through
        # without the caller needing to strip the index first. Symlinked under the
        # input's own basename (not FileA/FileB, whose localized names may differ)
        # so a ".tbi" pairs correctly with its bgzip file.
        ln -s ~{FileA} ~{basename(FileA)}
        ln -s ~{FileB} ~{basename(FileB)}
        ~{"ln -s " + FileAIndex + " " + basename(FileA) + ".tbi"}
        ~{"ln -s " + FileBIndex + " " + basename(FileB) + ".tbi"}

        bedtools intersect -a ~{basename(FileA)} -b ~{basename(FileB)} ~{Params} > ~{Prefix}.intersect.bed
    >>>

    runtime {
        docker: DockerImage
        memory: "8G"
        cpu: 2
        disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
        preemptible: 3
        maxRetries: 2
    }

    output {
        File Output = "~{Prefix}.intersect.bed"
    }
}
