#!/usr/bin/env python3
# Replaces the extracted-sequence ID (<indiv>-...asm_h1...minimap2_SVA_0001)
# with indiv, hap and sva_id columns; every other column is passed through.
import sys
import os
import polars as pl
c = pl.col

file = sys.argv[1]
delim = os.environ["DELIMITER"]

# Read as strings so passed-through values are written back byte-identical.
df = pl.read_csv(file, separator=delim, infer_schema=False)

df = df.select(
    c("ID").str.split("-").list.get(0).alias("indiv"),
    c("ID").str.extract("asm_(h[12])", 1).alias("hap"),
    c("ID").str.extract(r"minimap2_(SVA_\d{4})", 1).alias("sva_id"),
    pl.all().exclude("ID"),
)

df.write_csv(sys.stdout, separator=delim)
