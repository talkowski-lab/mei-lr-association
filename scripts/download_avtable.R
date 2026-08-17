library(dplyr)
library(readr)
library(AnVIL)
library(argparser)


## Downloads an entire Terra data table once, so that per-row lookups
## downstream (e.g. in a scatter) can read a local file instead of
## re-querying Terra's AVTable API once per row.

argv <- arg_parser("Download a full Terra data table to a local TSV") %>%
  add_argument("--avtable-name", help = "Terra data table (entity type) name, e.g. eQTL_COMB_9k") %>%
  add_argument("--workspace-namespace", help = "Terra workspace namespace (billing project) hosting the data table") %>%
  add_argument("--workspace-name", help = "Terra workspace name hosting the data table") %>%
  parse_args()

av_table <- avtable(argv$avtable_name, namespace = argv$workspace_namespace, name = argv$workspace_name)

write_tsv(av_table, "avtable.tsv")
