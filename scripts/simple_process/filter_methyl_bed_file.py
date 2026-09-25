#!/usr/bin/env python3
import sys
import os
import polars as pl
c = pl.col

file = sys.argv[1]
df = pl.read_csv(file, separator=os.environ["DELIMITER"], has_header=False)
df.columns = [
    "chrom", "5mC_start", "5mC_end", "5mC_score", "haplotype", "coverage",
    "est_5mC_site_count", "est_unmod_site_count", "discrete_5mC_score",
    "_chrom", "_start", "_end", "var_id", "strand", "var_start", "var_end"
]

df = df.with_columns(
  
    # Signed distance from the (pre-window-expansion) variant to the 5mC site,
    # oriented by strand so downstream is positive. dist1/dist2 measure from
    # each end of the variant to each end of the site; opposite signs mean the
    # site falls inside the variant, i.e. distance 0.
    dist1 = pl.when(c("strand") == "-").then(c("var_start") - c("5mC_end")).otherwise(c("5mC_start") - c("var_end")),
    dist2 = pl.when(c("strand") == "-").then(c("var_end") - c("5mC_start")).otherwise(c("5mC_end") - c("var_start"))
).with_columns(
    distance = pl.when(c("dist1").sign() * c("dist2").sign() == -1).then(0).otherwise(c("dist1").sign() * pl.min_horizontal(c("dist1").abs(), c("dist2").abs()))
  )

df = df.drop("dist1", "dist2", "_chrom", "_start", "_end", "var_start", "var_end").filter(c("distance").abs() <= 100)

df.write_csv(sys.stdout, separator=os.environ["DELIMITER"])


