#!/usr/bin/env python3
import sys
import os
import polars as pl
c = pl.col

file = sys.argv[1]
df = pl.read_csv(file, separator=os.environ["DELIMITER"])

df = df.select(
    c("ID").str.split("-").list.get(0).alias("indiv"),
    c("ID").str.extract("asm_(h[12])", 1).alias("hap"),
    c("ID").str.extract(r"minimap2_(SVA_\d{4})", 1).alias("sva_id"),
    "length_hexamer",
    "length_VNTR"
)

df.write_csv(sys.stdout, separator=os.environ["DELIMITER"])


