#!/usr/bin/env Rscript
## Render an Rmarkdown file to a report.
##
## The Rmd is expected to save any figures/tables it generates into
## --figures-dir / --tables-dir (paths relative to the working directory this
## script is run from) -- this script only creates those directories and
## renders the report; it does not inspect or move chunk output itself.

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
  library(jsonlite)
})

FORMAT_MAP <- list(
  pdf = "pdf_document",
  html = "html_document",
  word = "word_document",
  docx = "word_document"
)

option_list <- list(
  make_option("--rmd", type = "character",
              help = "Path to the .Rmd file to render."),
  make_option("--format", type = "character", default = "pdf",
              help = paste0("Output format: ", paste(names(FORMAT_MAP), collapse = ", "),
                             ". [default %default]")),
  make_option("--output-file", type = "character", default = NULL,
              help = "Rendered report filename. Defaults to 'report.<ext>'."),
  make_option("--figures-dir", type = "character", default = "figures",
              help = "Directory the Rmd writes figures to; created if missing. [default %default]"),
  make_option("--tables-dir", type = "character", default = "tables",
              help = "Directory the Rmd writes tables to; created if missing. [default %default]"),
  make_option("--params-json", type = "character", default = NULL,
              help = "Optional JSON file of named params passed to rmarkdown::render(params=)."),
  make_option("--continue-on-chunk-error", action = "store_true", default = FALSE,
              help = paste("Keep rendering past a chunk error instead of aborting. The error is shown",
                           "inline in the report (knitr's default behavior) and also appended to --error-log.")),
  make_option("--error-log", type = "character", default = "render_errors.log",
              help = "Path to append chunk error messages to when --continue-on-chunk-error is set. [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rmd)) stop("--rmd is required")

fmt <- FORMAT_MAP[[tolower(opt$format)]]
if (is.null(fmt)) {
  stop(sprintf("Unsupported --format '%s'; expected one of: %s",
               opt$format, paste(names(FORMAT_MAP), collapse = ", ")))
}
ext <- switch(fmt, pdf_document = "pdf", html_document = "html", word_document = "docx")
output_file <- if (is.null(opt[["output-file"]])) paste0("report.", ext) else opt[["output-file"]]

dir.create(opt[["figures-dir"]], recursive = TRUE, showWarnings = FALSE)
dir.create(opt[["tables-dir"]], recursive = TRUE, showWarnings = FALSE)

params_list <- NULL
if (!is.null(opt[["params-json"]])) {
  params_list <- jsonlite::fromJSON(opt[["params-json"]])
  if (length(params_list) == 0) params_list <- NULL
}

# Passing a plain format *name* to rmarkdown::render() means it always
# resets opts_chunk$error to FALSE part-way through its own setup, so
# error=TRUE must be attached to the output format *object* instead --
# rmarkdown applies output_format$knitr on top of that reset.
output_format <- switch(fmt,
  pdf_document = rmarkdown::pdf_document(),
  html_document = rmarkdown::html_document(),
  word_document = rmarkdown::word_document()
)

if (opt[["continue-on-chunk-error"]]) {
  error_log <- opt[["error-log"]]
  output_format$knitr$opts_chunk$error <- TRUE
  output_format$knitr$knit_hooks$error <- function(x, options) {
    label <- if (is.null(options$label)) "unnamed" else options$label
    msg <- paste(x, collapse = "\n")
    cat(sprintf("[%s] chunk '%s':\n%s\n\n", format(Sys.time()), label, msg),
        file = error_log, append = TRUE)
    sprintf("\n\n**Error in chunk `%s`:**\n\n```\n%s\n```\n\n", label, msg)
  }
}

rmarkdown::render(
  input = opt$rmd,
  output_format = output_format,
  output_file = output_file,
  output_dir = ".",
  knit_root_dir = normalizePath("."),
  params = params_list,
  envir = new.env()
)
