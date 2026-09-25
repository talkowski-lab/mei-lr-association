#!/usr/bin/env python3
import sys
import os
import polars as pl
c = pl.col

selected_svas = [
    "SVA_1049", # Disease relevant SVAs
    "SVA_0481", "SVA_1885", "SVA_1056", "SVA_0438", "SVA_1404", "SVA_1049", "SVA_2107", "SVA_1700", "SVA_1447", "SVA_1433", # Longest Avg Hex
    "SVA_1885", "SVA_0481", "SVA_0110", "SVA_0438", "SVA_1056", "SVA_1049", "SVA_0184", "SVA_1156", "SVA_1190", "SVA_2051", # Most variable Hex
    "SVA_1403", "SVA_0861", "SVA_0923", "SVA_1258", "SVA_1852", "SVA_0311", "SVA_0859", "SVA_0202", "SVA_0870", "SVA_1636", # Longest Avg VNTRs
    "SVA_0861", "SVA_1403", "SVA_0228", "SVA_0870", "SVA_0933", "SVA_0528", "SVA_1127", "SVA_1676", "SVA_1258", "SVA_1852"  # Most variable VNTR
]

file = sys.argv[1]
df = pl.read_csv(file, separator=os.environ["DELIMITER"])

df = df.select(
    c("ID").str.split("-").list.get(0).alias("indiv"),
    c("ID").str.extract("asm_(h[12])", 1).alias("hap"),
    c("ID").str.extract(r"minimap2_(SVA_\d{4})", 1).alias("sva_id"),
    "length_hexamer",
    "length_VNTR"
).filter(c("sva_id").is_in(selected_svas))

df.write_csv(sys.stdout, separator=os.environ["DELIMITER"])


