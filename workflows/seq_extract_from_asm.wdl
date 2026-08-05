version 1.0

workflow ExtractSeqFromASM {
  input {
    File hap1_bam
    File hap1_bam_index
    File hap2_bam
    File hap2_bam_index
    File region_file
    String Prefix
    String ImageTag = "latest"
  }


  call ExtractSeq as ExtractSeq1 {
    input: region_file = region_file, input_file = hap1_bam, input_file_index=hap1_bam_index, ImageTag = ImageTag
  }

  call ExtractSeq as ExtractSeq2 {
    input: region_file = region_file, input_file = hap2_bam, input_file_index=hap2_bam_index, ImageTag = ImageTag

  }

  call CombineSeqs {
    input: fa1 = ExtractSeq1.sva_seq, fa2 = ExtractSeq2.sva_seq, prefix = Prefix
  }


  output {
    File combined_sva_seqs = CombineSeqs.combined_sva_seqs
  }
}

task ExtractSeq {
  input {
    File region_file
    File input_file
    File input_file_index
    String ImageTag = "latest"
    Int? DiskGB
  }

  String name = basename(input_file, ".bam")
  String outfile = basename(input_file, ".bam") + "_sva_seq.fa"

  Int auto_disk_size = ceil(size(input_file, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    ln -s ~{input_file} input.bam
    ln -s ~{input_file_index} input.bam.bai

    fwd_file="~{name}_fwd.fasta"
    rev_file="~{name}_rev.fasta"

    awk '$5 == "+"' ~{region_file} | cut -f 1-3 > fwd_regions.bed
    awk '$5 == "+"' ~{region_file} | cut -f 7 > fwd_ids.txt
    awk '$5 == "-"' ~{region_file} | cut -f 1-3 > rev_regions.bed
    awk '$5 == "-"' ~{region_file} | cut -f 7 > rev_ids.txt

    # One samtools consensus call per strand (instead of one per region) via
    # --regions-file (samtools >= 1.24); each output record's header is renamed
    # positionally using the matching id list, since --regions-file itself has
    # no way to assign a custom per-region sequence name.
    rename_headers() {
      ids_file=$1
      awk -v ids_file="$ids_file" -v prefix="~{name}_" '
        BEGIN { i = 0; while ((getline id < ids_file) > 0) ids[++i] = id }
        /^>/ { n++; print ">" prefix ids[n]; next }
        { print }
      '
    }

    samtools consensus --regions-file fwd_regions.bed -m simple input.bam \
      | rename_headers fwd_ids.txt > $fwd_file

    samtools consensus --regions-file rev_regions.bed -m simple input.bam \
      | rename_headers rev_ids.txt > $rev_file

    python3 /scripts/combine_fastas.py $fwd_file <(seqtk seq -r $rev_file) > ~{outfile}

    if [ ! -s ~{outfile} ]; then
      echo "Error: No sequenes were extracted" >&2
      exit 1
    fi

  >>>

  output {
    File sva_seq = outfile
  }

  runtime {
    docker: "ayenkin1871/mei-lr-association-bioinformatics:" + ImageTag
    cpu: 4
    memory: "30 GiB"
    disks: "local-disk " + select_first([DiskGB, auto_disk_size]) + " HDD"
    bootDiskSizeGb: 10
    preemptible: 3
    maxRetries: 2
  }
}

task CombineSeqs {
  input {
    File fa1
    File fa2
    String prefix
   
  }
  String outfile = prefix + "_combined_sva_seqs.fa"
  
  command <<<
    cat ~{fa1} ~{fa2} > ~{outfile}
  >>>

  output {
    File combined_sva_seqs = outfile
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 2
    memory: "4 GiB"
    disks: "local-disk 10 HDD"
    bootDiskSizeGb: 10
    preemptible: 3
    maxRetries: 2
  }
}

