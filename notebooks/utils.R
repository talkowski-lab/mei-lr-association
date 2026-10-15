library(glue)

set_gcloud_bucket <- function(bucket) options(gcloud_bucket = bucket)


gs_download <- function(file, folder=".", force=FALSE) {
  if (!startsWith(file, "gs://")) {
    bucket <- getOption("gcloud_bucket")
    if (is.null(bucket)) {
      stop("Must provide gcloud url or default bucket")
    } else {
      file <- paste0("gs://", bucket, gsub("^/", "", file))
    }
  }
  basef <- basename(file)
  outf <- glue("{folder}/{basef}")
  if (file.exists(outf) && !force) {
    cat(glue("{outf} already exists, set force=TRUE to redownload"))
    return(outf)
  }
  command <- glue("gcloud storage cp {file} {folder}/{basef}")
  system(command = command)
  return(glue("{folder}/{basef}"))
}

gs_upload <- function(file, folder) {
  if (!startsWith(file, "gs://")) {
    bucket <- getOption("gcloud_bucket")
    if (is.null(bucket)) {
      stop("Must provide gcloud url or default bucket")
    } else {
      file <- paste0("gs://", bucket, gsub("^/", "", file))
    }
  }
  basef <- basename(file)
  if (endsWith(folder, "/")) {
    command <- glue("gcloud storage cp {file} {folder}{basef}")
  } else {
    command <- glue("gcloud storage cp {file} {folder}/{basef}")
  }
  system(command = command)
  return(glue("{folder}/{basef}"))
}


within_ranges <- function(v, range_starts, N) {
  sum(Reduce(`|`, lapply(range_starts, function(r) { v >= r & v < r + N}), init=FALSE))
}

n_intervals_cover <- function(vec, N, frac=0.95) {
  # Remove NA values
  vec <- vec[!is.na(vec)]
  threshold <- ceiling(frac * length(vec))
  unique_vals <- unique(vec)
  
  count <- 0
  final_intervals <- c()
  while (count < threshold) {
    curr_interval <- -N
    
    for (val in unique_vals) {
      wr_count <- within_ranges(vec, c(final_intervals, val), N)
      if (wr_count > count) {
        curr_interval <- val
        count <- wr_count
      }
    }
    final_intervals <- c(final_intervals, curr_interval)
    unique_vals <- unique_vals[unique_vals != curr_interval]
  }
  return(length(final_intervals))
}


download_and_read_batched_delim <- function(avtable, column_name, outfolder=".", n_batch=200,  transform_func_indiv=NULL, transform_func_batch=NULL, prefix_col=NULL, total_n=Inf, ...) {
  #TODO: Get rid of prefix col the files should have unique names
  
  if (!dir.exists(outfolder)) {
    dir.create(outfolder, recursive=TRUE)
  }
  
  n_batches <- ceiling(min(nrow(avtable), total_n) / n_batch)
  
  total_df <- lapply(0:(n_batches-1), function(i) {
    cat(glue("\r{i+1} out of {n_batches}"))
    batch_start <- i * n_batch + 1
    batch_end <- min(nrow(table), (i+1) * n_batch)
    
    files <- avtable[[column_name]][batch_start:batch_end]
    files <- files[!is.na(files)]
    if (length(files) == 0) {
      return(NULL)
    }
    
    if (!is.null(prefix_col)) {
      prefixes <- avtable[[prefix_col]][batch_start:batch_end][!is.na(files)]
      outfiles <- sapply(files, function(f) glue("{outfolder}/{prefixes}_{basename(f)}"), USE.NAMES=FALSE)
    } else {
      outfiles <- sapply(files, function(f) glue("{outfolder}/{basename(f)}"), USE.NAMES=FALSE)
    }
    
    download_files <- files[!file.exists(outfiles)]
    
    if (length(download_files) > 0) {
      command <- paste("gcloud storage cp ", paste(download_files, collapse=" "), outfolder)
      system(command)
    }
    print(outfiles)
    df <- lapply(outfiles, function(f) {
      tab <- read_delim(f, progress=FALSE, show_col_types=FALSE, ...)
      if (!is.null(transform_func_indiv)) {
        tab <- transform_func_indiv(tab, f)
      }
      return(tab)
    }) %>%
      bind_rows()
    
    if (!is.null(transform_func_batch)) {
      df <- transform_func(df)
    }
    return(df)
  }) %>%
    bind_rows()
  print(total_df)
  return(total_df)
}



read_one_parquet <- function(path) {
  message("Reading: ", path)
  
  out <- tryCatch(
    arrow::read_parquet(path),
    error = function(e) {
      message("Skipping file due to read error: ", path)
      return(NULL)
    }
  )
  
  if (is.null(out)) {
    return(NULL)
  }
  
  out <- as.data.frame(out)
  
  if (nrow(out) == 0) {
    message("File is empty: ", path)
    return(NULL)
  }
  
  out
}
