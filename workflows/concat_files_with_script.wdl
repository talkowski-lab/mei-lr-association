version 1.0

import "utils/concat_files.wdl" as ConcatFiles

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
## files, one giant ConcatenateDelim call pays for all of that localization
## serially on a single VM. If BatchSize is set, InputFiles is instead chunked
## into groups of that size, each chunk is localized/concatenated by its own
## ConcatenateDelim call (so localization is spread across many concurrent
## tasks instead of one), and the resulting per-batch files are concatenated
## again by a final ConcatenateDelim call (which is where GzipOutput is
## actually applied -- per-batch outputs are always left ungzipped so this
## final call can still parse them as text). If BatchSize is unset, this
## collapses back to the original single-call behavior with no scatter at all.

workflow ConcatenateAndProcessFiles {
  input {
    Array[File] InputFiles
    Boolean? HasHeader
    Boolean? GzipOutput
    File? RemoteScript
    String? EmbeddedScriptPath
    Boolean? AddIdColumn
    Array[String]? Ids
    String? NewIdColumnName
    String? Delimiter
    String? OutputName
    String? ImageTag 
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

      call ConcatFiles.ConcatenateDelim as ConcatenateBatch {
        input:
          InputFiles = batch_file,
          HasHeader = HasHeader,
          GzipOutput = false,
          RemoteScript = RemoteScript,
          EmbeddedScriptPath = EmbeddedScriptPath,
          AddIdColumn = AddIdColumn,
          Ids = batch_id_elem,
          NewIdColumnName = NewIdColumnName,
          Delimiter = Delimiter,
          OutputName = "batch_" + batch_idx + ".txt",
          ImageTag = ImageTag
      }
    }

    # Batch outputs already have (at most) one header each and their rows are
    # already ID-tagged, so no script/ID-column options are passed here --
    # this call only re-applies header de-dup across batches and GzipOutput.
    call ConcatFiles.ConcatenateDelim as MergeBatches {
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
    call ConcatFiles.ConcatenateDelim as ConcatenateSingle {
      input:
        InputFiles = InputFiles,
        HasHeader = HasHeader,
        GzipOutput = GzipOutput,
        RemoteScript = RemoteScript,
        EmbeddedScriptPath = EmbeddedScriptPath,
        AddIdColumn = AddIdColumn,
        Ids = Ids,
        NewIdColumnName = NewIdColumnName,
        Delimiter = Delimiter,
        OutputName = OutputName,
        ImageTag = ImageTag
    }
  }

  output {
    File ConcatenatedFile = select_first([MergeBatches.ConcatenatedFile, ConcatenateSingle.ConcatenatedFile])
  }
}

