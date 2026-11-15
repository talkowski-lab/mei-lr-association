version 1.0

## Lists the data tables in a Terra workspace using terra-notebook-utils.

workflow ListTerraTables {
    input {
        String Workspace
        String WorkspaceNamespace
        String ImageTag = "latest"
    }

    call ListTables {
        input:
            Workspace = Workspace,
            WorkspaceNamespace = WorkspaceNamespace,
            ImageTag = ImageTag
    }

    output {
        File TableList = ListTables.TableList
    }
}

task ListTables {
    input {
        String Workspace
        String WorkspaceNamespace
        String ImageTag = "latest"
        Int MemoryGB = 2
    }

    command <<<
        set -euo pipefail
        tnu table list \
            --workspace "~{Workspace}" \
            --workspace-namespace "~{WorkspaceNamespace}" \
            | tee tables.txt
    >>>

    output {
        File TableList = "tables.txt"
    }

    runtime {
        docker: "ayenkin1871/mei-lr-association-terra:" + ImageTag
        memory: "~{MemoryGB} GB"
        cpu: 1
        disks: "local-disk 10 HDD"
        preemptible: 2
    }
}
