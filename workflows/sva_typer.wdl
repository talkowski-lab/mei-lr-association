version 1.0

workflow sva_extract_workflow {
  input {
    File sva_fasta
    String docker_file = "ayenkin1871/sva_typer:latest"
  }

  call SVA_typer {
    input:  sva_seq_file = sva_fasta, docker_file = docker_file 
  }

  output {
    File repeat_positions = SVA_typer.repeat_positions
    File repeat_lengths = SVA_typer.repeat_lengths
  }
}


task SVA_typer {
  input {
    File sva_seq_file
    String docker_file
  }
  
  String base = basename(sva_seq_file, ".fa")
  String outpos = base + "_repeat_positions.txt"
  String outfile = base + "_repeat_lengths.txt"
  
  command <<<
    sva_typer --write-query-seq-state ~{sva_seq_file} > ~{outpos}
    python /usr/src/sva_typer/scripts/process_output.py ~{outpos} > ~{outfile}
  >>>

  output {
    File repeat_lengths = outfile
    File repeat_positions = outpos
  }

  runtime {
    docker: docker_file
    cpu: 4
    memory: "30 GiB"
    disks: "local-disk 50 HDD"
    bootDiskSizeGb: 10
    preemptible: 3
    maxRetries: 2
  }
}
