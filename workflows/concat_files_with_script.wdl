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
##
## BatchSize is optional and off by default. Cromwell localizes every File in
## InputFiles before a task's command even starts, so with thousands of small
## files, one giant ConcatenateFiles call pays for all of that localization
## serially on a single VM. If BatchSize is set, InputFiles is instead chunked
## into groups of that size, each chunk is localized/concatenated by its own
## ConcatenateFiles call (so localization is spread across many concurrent
## tasks instead of one), and the resulting per-batch files are concatenated
## again by a final ConcatenateFiles call (which is where GzipOutput is
## actually applied -- per-batch outputs are always left ungzipped so this
## final call can still parse them as text). If BatchSize is unset, this
## collapses back to the original single-call behavior with no scatter at all.

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
    Int? BatchSize
  }

  Int total_files = length(InputFiles)

  if (defined(BatchSize)) {
    Int effective_batch_size = select_first([BatchSize])
    Int n_batches = (total_files + effective_batch_size - 1) / effective_batch_size

    scatter (batch_idx in range(n_batches)) {
      # No slicing operator in WDL 1.0, so InputFiles[batch_idx*size : ...] isn't
      # available -- gather this batch's indices via select_all over a masked
      # range, then gather the actual Files (and Ids, if given) by those indices.
      # None of this touches file *contents*, so nothing gets localized here --
      # only the ConcatenateBatch call below localizes anything, and only the
      # ~effective_batch_size files that belong to its one batch.
      scatter (file_pos in range(total_files)) {
        # WDL 1.0 has no `None` literal, so an if-block (rather than a
        # then/else ternary) is the portable way to produce an Int? that's
        # only set when file_pos falls in this batch.
        if (file_pos >= batch_idx * effective_batch_size && file_pos < (batch_idx + 1) * effective_batch_size) {
          Int keep_pos = file_pos
        }
      }
      Array[Int] batch_indices = select_all(keep_pos)

      scatter (idx in batch_indices) {
        File batch_file = InputFiles[idx]
      }

      if (defined(Ids)) {
        scatter (idx in batch_indices) {
          String batch_id_elem = select_first([Ids])[idx]
        }
      }

      call ConcatenateFiles as ConcatenateBatch {
        input:
          InputFiles = batch_file,
          HasHeader = HasHeader,
          GzipOutput = false,
          RemoteScript = RemoteScript,
          EmbeddedScriptPath = EmbeddedScriptPath,
          AddIdColumn = AddIdColumn,
          Ids = batch_id_elem,
          IdColumnName = IdColumnName,
          Delimiter = Delimiter,
          OutputName = "batch_" + batch_idx + ".txt",
          ImageTag = ImageTag
      }
    }

    # Batch outputs already have (at most) one header each and their rows are
    # already ID-tagged, so no script/ID-column options are passed here --
    # this call only re-applies header de-dup across batches and GzipOutput.
    call ConcatenateFiles as MergeBatches {
      input:
        InputFiles = ConcatenateBatch.ConcatenatedFile,
        HasHeader = HasHeader,
        GzipOutput = GzipOutput,
        Delimiter = Delimiter,
        OutputName = OutputName,
        ImageTag = ImageTag
    }
  }

  if (!defined(BatchSize)) {
    call ConcatenateFiles as ConcatenateSingle {
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
  }

  output {
    File ConcatenatedFile = select_first([MergeBatches.ConcatenatedFile, ConcatenateSingle.ConcatenatedFile])
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

  # Ids is Array[String] (never File), so it's safe to write out here regardless
  # of localization timing. InputFiles is NOT handled this way -- see the manifest
  # built inline in the command block below, and the comment there for why.
  File id_manifest = write_lines(select_first([Ids, []]))

  # Only append .gz if OutputName doesn't already end with it, so callers who
  # already anticipated gzipping (e.g. OutputName = "concatenated.tsv.gz")
  # don't end up with a doubled ".gz.gz".
  Boolean output_already_has_gz_suffix = sub(OutputName, "\\.gz$", "") != OutputName
  String FinalOutputName = if GzipOutput && !output_already_has_gz_suffix then OutputName + ".gz" else OutputName
  Boolean AppendedGzSuffix = FinalOutputName != OutputName

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
    export OUT_PATH='~{FinalOutputName}'
    export ORIGINAL_OUTPUT_NAME='~{OutputName}'
    export APPENDED_GZ_SUFFIX='~{AppendedGzSuffix}'

    # Built inline (rather than via a `File manifest = write_lines(InputFiles)`
    # declaration outside command) so each line is InputFiles' *localized* local
    # path: a File only gets its Cromwell-localized path substituted when
    # referenced directly inside the command block, not in declarations
    # evaluated before it -- same gotcha as the RemoteScript "+" idiom above.
    cat > input_files_manifest.txt <<'MANIFEST_EOF'
~{sep="\n" InputFiles}
MANIFEST_EOF

    python3 <<'CODE'
import csv
import gzip
import io
import os
import subprocess
import sys

has_header = os.environ["HAS_HEADER"] == "true"
gzip_output = os.environ["GZIP_OUTPUT"] == "true"
add_id_column = os.environ["ADD_ID_COLUMN"] == "true"
id_col_name = os.environ["ID_COLUMN_NAME"]
delim = os.environ["DELIMITER"]
remote_script_path = os.environ.get("REMOTE_SCRIPT_PATH") or None
embedded_script_path = os.environ.get("EMBEDDED_SCRIPT_PATH") or None
out_path = os.environ["OUT_PATH"]  # already has .gz appended (if any) by the WDL's FinalOutputName

if os.environ["APPENDED_GZ_SUFFIX"] == "true":
    print(
        f"[info] GzipOutput is set and OutputName ({os.environ['ORIGINAL_OUTPUT_NAME']!r}) "
        f"had no .gz suffix -- appended one; writing to {out_path!r} instead",
        file=sys.stderr,
    )

with open("input_files_manifest.txt") as f:
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

if script_path is not None:
    print(f"[debug] preprocessing script: {script_path!r} exists={os.path.isfile(script_path)}", file=sys.stderr)
else:
    print("[debug] no preprocessing script set; using file contents as-is", file=sys.stderr)

# compresslevel=6 to match the standalone gzip CLI's default -- gzip.open()'s
# own default (9, max compression) is meaningfully slower for little ratio gain.
out_f_cm = gzip.open(out_path, "wt", newline="", compresslevel=6) if gzip_output else open(out_path, "wt", newline="")
writer = None
with out_f_cm as out_f:
    for i, fp in enumerate(files):
        if script_path is not None:
            result = subprocess.run([script_path, fp], capture_output=True)
            if result.returncode != 0:
                print(f"[error] {script_path} failed on {fp!r} (exit {result.returncode})", file=sys.stderr)
                print(result.stderr.decode(errors="replace"), file=sys.stderr)
                result.check_returncode()
            text = result.stdout.decode()
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
