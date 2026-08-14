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
