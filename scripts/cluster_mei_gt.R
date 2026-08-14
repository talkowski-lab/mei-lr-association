library(data.table)
library(glue)
library(dplyr)
library(readr)
library(tidyr)
library(argparser)
library(stringr)
library(digest)


argv <- arg_parser("Cluster MEI GTs from a VCF") %>%
  add_argument("--mei-gt-mat", "The matrix containing GTs for MEIs") %>%
  add_argument("--mei-cluster-bed", "Bed file containing MEI information and cluster #") %>%
  add_argument("--prefix", "Output prefix") %>%
  parse_args()

infile <- argv$mei_gt_mat
outprefix <- "tidy_gt_"

# read only header
header <- names(fread(infile, nrows = 0))
header <- str_replace_all(header, c("^#"="", "\\[\\d+\\]"="", ":GT"=""))
indivs <- header[5:length(header)]

skip <- 1L     # start after header row
chunk_size <- 10000
chunk_i <- 0

repeat {
  chunk <- fread(
    infile,
    skip = skip,
    nrows = chunk_size,
    header = FALSE,
    col.names = header,
    colClasses = c(
      c("character", "integer", "character", "character"),
      rep("character", length(indivs))
    ),
    na.strings = c("./.", ".")
  )

  if (nrow(chunk) == 0) break
  
  chunk <- tibble(chunk)
  
  tidy_chunk <- chunk %>%
    pivot_longer(cols=-c(CHROM, POS, REF, ALT), 
                 names_to="indiv", 
                 values_to="GT", 
                 values_drop_na=TRUE) 
  
  chunk_filename <- paste0(outprefix, chunk_i, ".txt")
  tidy_chunk %>%
    write_delim(chunk_filename, delim="\t")
  

  skip <- skip + nrow(chunk)
  chunk_i <- chunk_i + 1
  rm(chunk, tidy_chunk)
  gc()
}

bed_df <- read_delim(argv$mei_cluster_bed) %>%
  mutate(ALT_md5 = sapply(ALT, digest, algo="md5", serialize=FALSE))

#Now generate final table, joining with length table
gt_tidy <- lapply(0:(chunk_i-1), function(i) {
  chunk_filename <- paste0(outprefix, i, ".txt")
  df <- read_delim(chunk_filename) %>%
    mutate(ALT_md5 = sapply(ALT, digest, algo="md5", serialize=FALSE)) %>% 
    select(CHROM, POS, ALT_md5, indiv, GT) 
}) %>%
  bind_rows() %>%
  left_join(bed_df %>% select(chrom, start, end, ALT_md5, cluster_id), by=c(CHROM="chrom", POS="start", "ALT_md5"))
  
## Check for duplicate genotypes across haplotypes
gt_tidy <- gt_tidy %>%
  mutate(
    hap1 = GT %in% c("1", "1/1", "1|1", "1|0"),
    hap2 = case_when(
      GT %in% c("1/1", "1|1", "0|1") ~ TRUE,
      GT %in% c("0/0", "1|0") ~ FALSE,
      .default=NA
  )
)


gt_dups <- gt_tidy %>%
  group_by(cluster_id, indiv) %>%
  filter(sum(hap1) > 1 | sum(hap2) > 1) %>%
  arrange(CHROM, cluster_id, indiv)

dup_clusters <- gt_dups$cluster_id %>% unique()

gt_tidy %>% 
  filter(!(cluster_id %in% dup_clusters)) %>%
  write_tsv(glue("{argv$prefix}_gt_tidy_clustered.txt"))

gt_tidy %>% 
  filter(cluster_id %in% dup_clusters) %>%
  write_tsv(glue("{argv$prefix}_gt_tidy_dupclusters.txt"))
