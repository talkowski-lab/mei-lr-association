version 1.0

## Downloads an entire Terra data table once, so callers doing many per-row
## lookups (e.g. in a scatter) can read a local file instead of each hitting
## Terra's AVTable API separately.
task DownloadAVTable {
  input {
    String AVTableName
    String WorkspaceNamespace
    String WorkspaceName
    String ImageTag = "latest"
    Int MemoryGB = 4
    Int? DiskGB
  }

  command <<<
    set -euo pipefail

    Rscript /scripts/download_avtable.R \
        --avtable-name ~{AVTableName} \
        --workspace-namespace ~{WorkspaceNamespace} \
        --workspace-name ~{WorkspaceName}
  >>>

  runtime {
    docker: "ayenkin1871/mei-lr-association-r_analysis:" + ImageTag
    memory: MemoryGB + " GB"
    cpu: 2
    disks: "local-disk " + select_first([DiskGB, 20]) + " SSD"
    preemptible: 3
    maxRetries: 2
  }

  output {
    File AVTableData = "avtable.tsv"
  }
}
