version 1.0

## Exports every data table of each workspace in Workspaces
## ("<namespace>/<workspace>") to "<namespace>.<workspace>.<table>.tsv".
## Runs with the task's own credentials (the pet service account on Terra),
## which need read access to each workspace.
task ExportTerraTables {
  input {
    Array[String] Workspaces
    String ImageTag = "latest"
    Int MemoryGB = 4
    Int DiskGB = 20
  }

  meta {
    volatile: true
  }

  command <<<
    set -euo pipefail
    python3 /scripts/export_terra_tables.py \
        --out-dir terra_tables \
        ~{sep=" " prefix("--workspace ", Workspaces)}
  >>>

  runtime {
    docker: "ayenkin1871/mei-lr-association-terra:" + ImageTag
    memory: MemoryGB + " GB"
    cpu: 1
    disks: "local-disk " + DiskGB + " HDD"
    preemptible: 2
    maxRetries: 1
  }

  output {
    Array[File] Tables = glob("terra_tables/*.tsv")
  }
}
