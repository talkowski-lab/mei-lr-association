version 1.0

import "utils/split_files.wdl" as SplitFiles
import "utils/concat_files.wdl" as ConcatFiles

workflow sva_extract_workflow {
  input {
    File SVAFasta
    String ImageTag = "latest"
    Int? BatchSize
  }
  String base = basename(SVAFasta, ".fa")

  if (defined(BatchSize)) {
    Int effective_batch_size = select_first([BatchSize])

    call SplitFiles.SplitFasta as SplitSVASeqs {
      input: fasta = SVAFasta, records_per_file = effective_batch_size
    }

    scatter (fasta in SplitSVASeqs.split_fastas) {

      String batch_prefix = basename(fasta, ".fa")
      call SVA_typer as SVA_typer_Batch {
        input: sva_seq_file = fasta, ImageTag = ImageTag, Prefix = batch_prefix
      }
    }
    String outpos = base + "_repeat_positions.txt"
    String outfile = base + "_repeat_lengths.txt"

    call ConcatFiles.ConcatenateDelim as ConcatLengths {
      input:
        InputFiles = SVA_typer_Batch.repeat_lengths,
        OutputName = outfile
    }

    call ConcatFiles.ConcatenateDelim as ConcatPositions {
      input:
        InputFiles = SVA_typer_Batch.repeat_positions,
        OutputName = outpos
    }

  }

  if (!defined(BatchSize)) {
    call SVA_typer as SVA_typer_Single {
      input:  sva_seq_file = SVAFasta, ImageTag = ImageTag, Prefix = base
    }

  }


  output {
    File repeat_positions = select_first([SVA_typer_Single.repeat_positions, ConcatPositions.ConcatenatedFile])
    File repeat_lengths = select_first([SVA_typer_Single.repeat_lengths, ConcatLengths.ConcatenatedFile])
  }
}


task SVA_typer {
  input {
    File sva_seq_file
    String Prefix
    String ImageTag = "latest"
    Int CPU = 2
    Int MemoryGB = 4
    Int? DiskGB
  }

  String outpos = Prefix + "_repeat_positions.txt"
  String outfile = Prefix + "_repeat_lengths.txt"

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
