version 1.0

workflow GzipFile {
    input {
        File InputFile
        String ImageTag = "latest"
    }

    call Gzip {
        input:
            InputFile = InputFile,
            ImageTag = ImageTag
    }

    output {
        File GzippedFile = Gzip.Output
    }
}

task Gzip {
    input {
        File InputFile
        String ImageTag = "latest"
    }

    Int disk_size = ceil(size(InputFile, "GB") * 2) + 10

    command <<<
        set -euo pipefail

        gzip -c ~{InputFile} > ~{basename(InputFile)}.gz
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-bioinformatics:" + ImageTag
        memory: "4G"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File Output = basename(InputFile) + ".gz"
    }
}
