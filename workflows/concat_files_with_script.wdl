version 1.0

## Concatenates a list of files vertically, optionally running each file
## through a preprocessing script first and/or tagging rows with a per-file
## ID column (this ID-column logic was merged in from, and supersedes, the
## now-removed ConcatenateColumnFiles workflow in concat_files.wdl).
##
## RemoteScript (a File, e.g. a gs:// path -- localized by Cromwell like any
## other File input) and EmbeddedScriptPath (a path already baked into the
## Docker image, e.g. "/scripts/foo.py") are mutually exclusive; if both are
## set, RemoteScript wins and EmbeddedScriptPath is ignored. If neither is
## set, each file's raw contents are used as-is. Whichever script runs is
## invoked once per file as `<script> <file>` and must write the transformed
## file's contents to stdout.
##
## Typical Terra wiring for the ID column:
##   InputFiles = this.samples.my_file_column   (Array[File])
##   Ids        = this.samples.sample_id         (Array[String], optional)
## If Ids is omitted while AddIdColumn is set, each row is tagged with the
## basename of its (post-script) source file instead.

workflow ConcatenateAndProcessFiles {
  input {
    Array[File] InputFiles
    Boolean HasHeader = true
    Boolean GzipOutput = false
    File? RemoteScript
    String? EmbeddedScriptPath
    Boolean AddIdColumn = false
    Array[String]? Ids
    String IdColumnName = "source_id"
    String Delimiter = "\t"
    String OutputName = "concatenated.tsv"
    String ImageTag = "latest"
  }

  call ConcatenateFiles {
    input:
      InputFiles = InputFiles,
      HasHeader = HasHeader,
      GzipOutput = GzipOutput,
      RemoteScript = RemoteScript,
      EmbeddedScriptPath = EmbeddedScriptPath,
      AddIdColumn = AddIdColumn,
      Ids = Ids,
      IdColumnName = IdColumnName,
      Delimiter = Delimiter,
      OutputName = OutputName,
      ImageTag = ImageTag
  }

  output {
    File ConcatenatedFile = ConcatenateFiles.ConcatenatedFile
  }
}

task ConcatenateFiles {
  input {
    Array[File] InputFiles
    Boolean HasHeader
    Boolean GzipOutput
    File? RemoteScript
    String? EmbeddedScriptPath
    Boolean AddIdColumn
    Array[String]? Ids
    String IdColumnName
    String Delimiter
    String OutputName
    String ImageTag = "latest"
    Int CPU = 2
    Int MemoryGB = 4
    Int? DiskGB
  }

  # Manifest files avoid any shell-quoting issues with long/odd file paths.
  File manifest = write_lines(InputFiles)
  File id_manifest = write_lines(select_first([Ids, []]))

  String FinalOutputName = if GzipOutput then OutputName + ".gz" else OutputName

  Int auto_disk_size = ceil((size(InputFiles, "GiB") + size(RemoteScript, "GiB")) * 3) + 10

  command <<<
    set -euo pipefail

    # RemoteScript is referenced only via this "+" placeholder idiom (never
    # assigned to a separate File/String declaration outside command) so
    # Cromwell substitutes its localized local path here; declaring it
    # earlier would keep the raw gs:// URI instead.
    export REMOTE_SCRIPT_PATH='~{"" + RemoteScript}'
    export EMBEDDED_SCRIPT_PATH='~{select_first([EmbeddedScriptPath, ""])}'
    export HAS_HEADER='~{HasHeader}'
    export GZIP_OUTPUT='~{GzipOutput}'
    export ADD_ID_COLUMN='~{AddIdColumn}'
    export ID_COLUMN_NAME='~{IdColumnName}'
    export DELIMITER='~{Delimiter}'
    export OUT_PATH='~{OutputName}'

    python3 <<'CODE'
import csv
import gzip
import io
import os
import subprocess

has_header = os.environ["HAS_HEADER"] == "true"
gzip_output = os.environ["GZIP_OUTPUT"] == "true"
add_id_column = os.environ["ADD_ID_COLUMN"] == "true"
id_col_name = os.environ["ID_COLUMN_NAME"]
delim = os.environ["DELIMITER"]
remote_script_path = os.environ.get("REMOTE_SCRIPT_PATH") or None
embedded_script_path = os.environ.get("EMBEDDED_SCRIPT_PATH") or None
out_path = os.environ["OUT_PATH"]

with open("~{manifest}") as f:
    files = [line.strip() for line in f if line.strip()]

with open("~{id_manifest}") as f:
    ids = [line.strip() for line in f if line.strip()]

if ids and len(ids) != len(files):
    raise ValueError(
        f"Length of Ids ({len(ids)}) does not match length of InputFiles ({len(files)})"
    )

# RemoteScript and EmbeddedScriptPath are mutually exclusive; RemoteScript
# wins if both are set. Neither set means no preprocessing at all.
script_path = None
if remote_script_path:
    script_path = remote_script_path
    os.chmod(script_path, 0o755)  # localized GCS objects don't carry the exec bit
elif embedded_script_path:
    script_path = embedded_script_path

open_fn = gzip.open if gzip_output else open
writer = None
with open_fn(out_path, "wt", newline="") as out_f:
    for i, fp in enumerate(files):
        if script_path is not None:
            text = subprocess.run([script_path, fp], check=True, capture_output=True).stdout.decode()
        else:
            with open(fp) as in_f:
                text = in_f.read()

        rows = list(csv.reader(io.StringIO(text), delimiter=delim))
        if not rows:
            continue

        if has_header:
            header, data_rows = rows[0], rows[1:]
        else:
            header, data_rows = None, rows

        row_id = ids[i] if ids else os.path.basename(fp)

        if writer is None:
            writer = csv.writer(out_f, delimiter=delim)
            if header is not None:
                out_header = [id_col_name] + header if add_id_column else header
                writer.writerow(out_header)

        for row in data_rows:
            writer.writerow([row_id] + row if add_id_column else row)
CODE
  >>>

  output {
    File ConcatenatedFile = "~{FinalOutputName}"
  }

  runtime {
    docker: "ayenkin1871/mei-lr-association-python_general:" + ImageTag
    cpu: CPU
    memory: MemoryGB + " GB"
    disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
    preemptible: 3
    maxRetries: 2
  }
}
