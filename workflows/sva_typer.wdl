version 1.0

workflow sva_extract_workflow {
  input {
    File sva_fasta
    String ImageTag = "latest"
  }

  call SVA_typer {
    input:  sva_seq_file = sva_fasta, ImageTag = ImageTag
  }

  output {
    File repeat_positions = SVA_typer.repeat_positions
    File repeat_lengths = SVA_typer.repeat_lengths
  }
}


task SVA_typer {
  input {
    File sva_seq_file
    String ImageTag = "latest"
    Int CPU = 2
    Int MemoryGB = 4
    Int? DiskGB
  }

  String base = basename(sva_seq_file, ".fa")
  String outpos = base + "_repeat_positions.txt"
  String outfile = base + "_repeat_lengths.txt"

  Int auto_disk_size = ceil(size(sva_seq_file, "GB") * 2) + 10

  command <<<
    sva_typer --write-query-seq-state --hmm-behavior-n use ~{sva_seq_file} > ~{outpos}
    python /usr/src/sva_typer/scripts/process_output.py ~{outpos} > ~{outfile}
  >>>

  output {
    File repeat_lengths = outfile
    File repeat_positions = outpos
  }

  runtime {
    docker: "ayenkin1871/sva_typer:" + ImageTag
    cpu: CPU
    memory: MemoryGB + " GB"
    disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " SSD"
    bootDiskSizeGb: 10
    preemptible: 3
    maxRetries: 2
  }
}
