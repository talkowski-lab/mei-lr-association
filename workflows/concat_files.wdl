version 1.0

## Concatenates a set of delimited text files (e.g. a Terra data-table column
## of Array[File]) vertically into a single file, optionally prepending an
## ID column that records which source file (or entity) each row came from.
##
## Typical Terra wiring:
##   input_files = this.samples.my_file_column   (Array[File])
##   ids         = this.samples.sample_id         (Array[String], optional)
##
## If `ids` is omitted, each row is tagged with the basename of its source
## file instead.

workflow ConcatenateColumnFiles {
  input {
    Array[File] input_files          # the column of files to stack
    Array[String]? ids               # optional, same length/order as input_files
    Boolean add_id_column = true     # set false to just concatenate, no ID column
    String id_column_name = "source_id"
    Boolean has_header = true        # true if every input file has a header row
    String delimiter = "\t"          # "\t" for TSV, "," for CSV
    String output_name = "concatenated.tsv"
  }

  call ConcatenateFiles {
    input:
      input_files    = input_files,
      ids            = ids,
      add_id_column  = add_id_column,
      id_column_name = id_column_name,
      has_header     = has_header,
      delimiter      = delimiter,
      output_name    = output_name
  }

  output {
    File concatenated_file = ConcatenateFiles.concatenated_file
  }
}

task ConcatenateFiles {
  input {
    Array[File] input_files
    Array[String]? ids
    Boolean add_id_column
    String id_column_name
    Boolean has_header
    String delimiter
    String output_name

    Int? disk_size_gb
    Int mem_gb = 4
  }

  # Manifest files avoid any shell-quoting issues with long/odd file paths.
  File manifest = write_lines(input_files)
  File id_manifest = write_lines(select_first([ids, []]))

  Int auto_disk_size = ceil(size(input_files, "GiB") * 3) + 10

  command <<<
    set -euo pipefail

    python3 <<'CODE'
import csv
import os

with open("~{manifest}") as f:
    files = [line.strip() for line in f if line.strip()]

with open("~{id_manifest}") as f:
    ids = [line.strip() for line in f if line.strip()]

add_id = ~{if add_id_column then "True" else "False"}
has_header = ~{if has_header then "True" else "False"}
delim = "~{delimiter}"
id_col_name = "~{id_column_name}"
out_path = "~{output_name}"

if ids and len(ids) != len(files):
    raise ValueError(
        f"Length of ids ({len(ids)}) does not match length of input_files ({len(files)})"
    )

writer = None
with open(out_path, "w", newline="") as out_f:
    for i, fp in enumerate(files):
        row_id = ids[i] if ids else os.path.basename(fp)

        with open(fp, newline="") as in_f:
            reader = csv.reader(in_f, delimiter=delim)
            rows = list(reader)

        if not rows:
            continue

        if has_header:
            header, data_rows = rows[0], rows[1:]
        else:
            header, data_rows = None, rows

        if writer is None:
            writer = csv.writer(out_f, delimiter=delim)
            if header is not None:
                out_header = [id_col_name] + header if add_id else header
                writer.writerow(out_header)

        for row in data_rows:
            writer.writerow([row_id] + row if add_id else row)

print(f"Wrote {out_path}")
CODE
  >>>

  output {
    File concatenated_file = output_name
  }

  runtime {
    docker: "python:3.11-slim"
    memory: "~{mem_gb} GB"
    disks: "local-disk " + select_first([disk_size_gb, auto_disk_size]) + " HDD"
    cpu: 1
  }
}

