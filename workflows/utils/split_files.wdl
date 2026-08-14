version 1.0
## Splits a multi-record FASTA file into several smaller FASTA files,
## each containing up to `records_per_file` records.

task SplitFasta {
    input {
        File fasta
        Int records_per_file
        String output_prefix = "batch"
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
def read_fasta_records(path):
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, seq_lines

records_per_file = ~{records_per_file}
prefix = "~{output_prefix}"

file_index = 0
count = 0
out = None

for header, seq_lines in read_fasta_records("~{fasta}"):
    if count % records_per_file == 0:
        if out:
            out.close()
        file_index += 1
        out = open(f"{prefix}.{file_index}.fa", "w")
    out.write(header + "\n")
    if seq_lines:
        out.write("\n".join(seq_lines) + "\n")
    count += 1

if out:
    out.close()

print(f"Wrote {count} records to {file_index} file(s).")
CODE
    >>>

    output {
        Array[File] split_fastas = glob("~{output_prefix}.*.fa")
    }

    runtime {
        docker: "python:3.11-slim"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 20 SSD"
    }
}

## Splits a delimited file into several smaller files of up to `rows_per_file`
## rows each. If `has_header` is set, the first line is treated as a header
## and repeated in every output file (and excluded from the row count).
task SplitFileByRowCount {
    input {
        File input_file
        Int rows_per_file
        Boolean has_header = false
        String output_prefix = "batch"
    }

    command <<<
        set -euo pipefail
        shopt -s nullglob

        if [ "~{has_header}" = "true" ]; then
            header_line=$(head -n 1 ~{input_file})
            mkdir -p tmp_parts
            tail -n +2 ~{input_file} | split --numeric-suffixes=1 --suffix-length=6 -l ~{rows_per_file} - tmp_parts/part.
            for part in tmp_parts/part.*; do
                suffix=${part##*.}
                { echo "$header_line"; cat "$part"; } > ~{output_prefix}.${suffix}
            done
        else
            split --numeric-suffixes=1 --suffix-length=6 -l ~{rows_per_file} ~{input_file} ~{output_prefix}.
        fi
    >>>

    output {
        Array[File] split_files = glob("~{output_prefix}.*")
    }

    runtime {
        docker: "python:3.11-slim"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 20 SSD"
    }
}

## Splits a delimited file (with a header) into several smaller files, each
## containing up to `groups_per_file` distinct values of `group_column`. The
## input must already be sorted by that column, so each group's rows are
## contiguous -- a group boundary is detected purely by the column value
## changing from the previous row. The header is repeated in every output
## file.
task SplitFileByGroupCount {
    input {
        File input_file
        String group_column
        Int groups_per_file
        String delimiter = "\t"
        String output_prefix = "batch"
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import csv

delimiter = "~{delimiter}"
group_column = "~{group_column}"
groups_per_file = ~{groups_per_file}
prefix = "~{output_prefix}"

with open("~{input_file}", newline="") as f:
    reader = csv.reader(f, delimiter=delimiter)
    header = next(reader)
    if group_column not in header:
        raise SystemExit(f"Column {group_column!r} not found in header: {header}")
    group_idx = header.index(group_column)

    file_index = 0
    groups_seen = 0
    last_group = None
    out = None
    writer = None

    for row in reader:
        group = row[group_idx]
        if group != last_group:
            last_group = group
            groups_seen += 1
            if writer is None or groups_seen > groups_per_file:
                if out is not None:
                    out.close()
                file_index += 1
                out = open(f"{prefix}.{file_index}.txt", "w", newline="")
                writer = csv.writer(out, delimiter=delimiter)
                writer.writerow(header)
                groups_seen = 1
        writer.writerow(row)

    if out is not None:
        out.close()

print(f"Wrote {file_index} file(s).")
CODE
    >>>

    output {
        Array[File] split_files = glob("~{output_prefix}.*.txt")
    }

    runtime {
        docker: "python:3.11-slim"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 20 SSD"
    }
}
