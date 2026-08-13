version 1.0

workflow QueryVCF {
    input {
        File VCF
        File? VCFIndex
        String QueryFormat
        String? IncludeExpr
        String? ExcludeExpr
        String OutputName = "query_output.tsv"
        String ImageTag = "latest"
    }

    call BcftoolsQuery {
        input:
            VCF = VCF,
            VCFIndex = VCFIndex,
            QueryFormat = QueryFormat,
            IncludeExpr = IncludeExpr,
            ExcludeExpr = ExcludeExpr,
            OutputName = OutputName,
            ImageTag = ImageTag
    }

    output {
        File QueryOutput = BcftoolsQuery.Output
    }
}

task BcftoolsQuery {
    input {
        File VCF
        File? VCFIndex
        String QueryFormat
        String? IncludeExpr
        String? ExcludeExpr
        String OutputName
        String ImageTag = "latest"
    }

    command <<<
        set -euo pipefail

        bcftools query \
            ~{"-i '" + IncludeExpr + "'"} \
            ~{"-e '" + ExcludeExpr + "'"} \
            -f '~{QueryFormat}' \
            ~{VCF} > ~{OutputName}
    >>>

    runtime {
        docker: "ayenkin1871/mei-lr-association-bioinformatics:" + ImageTag
        memory: "8G"
        cpu: 2
        disks: "local-disk 50 HDD"
    }

    output {
        File Output = OutputName
    }
}
