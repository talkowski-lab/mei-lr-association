version 1.0

## Renders an Rmarkdown file into a report, plus whatever figures/tables the
## Rmd itself writes into FiguresDir/TablesDir. The Rmd is copied to the task's
## working directory (along with any InputData and RScripts files) before
## rendering, so it can refer to those files and to FiguresDir/TablesDir by
## relative path. Bucket (if given) is exported as the GCLOUD_BUCKET env var,
## for Rmds that build gs:// URIs from it (e.g. via gs_download/gs_upload
## helpers) instead of relying on an AnVIL/Terra workspace. If OutputDir (a
## gs:// folder) is given, the report, figures and tables are also copied to
## OutputDir/<Rmd basename without .Rmd>/{report.*,figures/,tables/}.

workflow RenderRmarkdownReport {
  input {
    File RmdFile
    Array[File] InputData = []
    Array[File] RScripts = []
    String FiguresDir = "figures"
    String TablesDir = "tables"
    String Format = "pdf"
    String DockerImage = "ayenkin1871/mei-lr-association-r_analysis:latest"
    Map[String, String] RenderParams = {}
    Boolean ContinueOnChunkError = false
    String? Bucket
    String? OutputDir
    Int MemoryGB = 8
    Int CPU = 2
    Int? DiskGB
  }

  call RenderReport {
    input:
      RmdFile = RmdFile,
      InputData = InputData,
      RScripts = RScripts,
      FiguresDir = FiguresDir,
      TablesDir = TablesDir,
      Format = Format,
      DockerImage = DockerImage,
      RenderParams = RenderParams,
      ContinueOnChunkError = ContinueOnChunkError,
      Bucket = Bucket,
      OutputDir = OutputDir,
      MemoryGB = MemoryGB,
      CPU = CPU,
      DiskGB = DiskGB
  }

  output {
    File Report = RenderReport.Report
    Array[File] Figures = RenderReport.Figures
    Array[File] Tables = RenderReport.Tables
    File? ErrorLog = RenderReport.ErrorLog
  }
}

task RenderReport {
  input {
    File RmdFile
    Array[File] InputData
    Array[File] RScripts
    String FiguresDir
    String TablesDir
    String Format
    String DockerImage
    Map[String, String] RenderParams
    Boolean ContinueOnChunkError
    String? Bucket
    String? OutputDir
    Int MemoryGB
    Int CPU
    Int? DiskGB
  }

  # Falls back to "out" (which the R script will reject) for an unrecognized format,
  # so the task fails with the script's clear error instead of a silent bad extension.
  String Ext = if Format == "pdf" then "pdf"
               else if Format == "html" then "html"
               else if Format == "word" || Format == "docx" then "docx"
               else "out"
  String OutputFile = "report." + Ext
  String NotebookName = basename(basename(RmdFile, ".Rmd"), ".rmd")

  File params_json = write_json(RenderParams)

  Int auto_disk_gb = ceil(size(RmdFile, "GB") + size(InputData, "GB") + size(RScripts, "GB")) * 2 + 20

  command <<<
    set -euo pipefail

    mkdir -p ~{FiguresDir} ~{TablesDir}

    export FIGURES_DIR=~{FiguresDir}
    export TABLES_DIR=~{TablesDir}
    export GCLOUD_BUCKET=~{select_first([Bucket, ""])}

    cp ~{RmdFile} report.Rmd

    # write_lines(Array[File]) can emit original (e.g. gs://) paths instead of the
    # localized on-VM paths, so the arrays are interpolated directly here instead --
    # that's the one substitution Cromwell always resolves to the localized path.
    for f in ~{sep=" " InputData}; do
        cp "$f" .
    done

    for f in ~{sep=" " RScripts}; do
        cp "$f" .
    done

    Rscript /scripts/render_rmarkdown.R \
        --rmd report.Rmd \
        --format ~{Format} \
        --output-file ~{OutputFile} \
        --figures-dir ~{FiguresDir} \
        --tables-dir ~{TablesDir} \
        --params-json ~{params_json} \
        ~{true="--continue-on-chunk-error" false="" ContinueOnChunkError}

    OUTPUT_DIR="~{select_first([OutputDir, ""])}"
    if [ -n "$OUTPUT_DIR" ]; then
        DEST="${OUTPUT_DIR%/}/~{NotebookName}"
        gcloud storage cp ~{OutputFile} "$DEST/"
        # Copies directory contents rather than the dir itself so a custom
        # FiguresDir/TablesDir still lands under figures/ and tables/; an empty
        # glob would make cp fail.
        if [ -n "$(ls -A ~{FiguresDir})" ]; then
            gcloud storage cp -r ~{FiguresDir}/* "$DEST/figures/"
        fi
        if [ -n "$(ls -A ~{TablesDir})" ]; then
            gcloud storage cp -r ~{TablesDir}/* "$DEST/tables/"
        fi
    fi
  >>>

  runtime {
    docker: DockerImage
    memory: MemoryGB + " GB"
    cpu: CPU
    disks: "local-disk " + select_first([DiskGB, auto_disk_gb]) + " SSD"
    preemptible: 2
    maxRetries: 1
  }

  output {
    File Report = OutputFile
    Array[File] Figures = glob(FiguresDir + "/*")
    Array[File] Tables = glob(TablesDir + "/*")
    File? ErrorLog = "render_errors.log"
  }
}
