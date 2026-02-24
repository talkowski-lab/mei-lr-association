import hail as hl
import subprocess
import pandas as pd

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--variant-list")
parser.add_argument("--sample-list")
parser.add_argument("--gs-vds-uri")
args = parser.parse_args()

hl.init()
hl.default_reference = "GRCh38"

with open(args.variant_list) as inf:
    variant_list = [line.strip() for line in inf]
    
with open(args.sample_list) as inf:
    user_ids = [line.strip() for line in inf]

intervals = []
for interval_str in variant_list:
    chrom, positions = interval_str.split(':')
    start = int(positions.split('_')[0])
    intervals.append({
        'interval': interval_str,
        'chrom': chrom,
        'start': start,
        'end': start + 1
    })


df = pd.DataFrame(intervals).sort_values(["chrom", "start"]).reset_index(drop=True)

binned_intervals = []
rows = list(df.itertuples(index=False))
current_chrom, current_start, current_end = rows[0].chrom, rows[0].start, rows[0].end

for row in rows[1:]:
    if row.chrom == current_chrom and row.end < current_start + 2_000_000:
        current_end = row.end
    else:
        binned_intervals.append(f"{current_chrom}:{current_start}-{current_end}")
        current_chrom, current_start, current_end = row.chrom, row.start, row.end

binned_intervals.append(f"{current_chrom}:{current_start}-{current_end}")

t = hl.import_table(args.variant_list, no_header=True, delimiter=r"\s+").rename({'f0': 'id'})
left_right = t.id.split(':')
rhs = left_right[1].split('_')
t = t.annotate(
    contig=left_right[0],
    pos=hl.int(rhs[0]),
    ref=rhs[1],
    alt=rhs[2],
)
t = t.annotate(
    locus=hl.locus(t.contig, t.pos, reference_genome="GRCh38"),
    alleles=[t.ref, t.alt]
).key_by('locus', 'alleles')

# Input parameters
vds = hl.vds.read_vds(args.gs_vds_uri)
vds_f = hl.vds.filter_intervals(
    vds,
    [hl.parse_locus_interval(x, reference_genome="GRCh38")
     for x in binned_intervals])

vds_f = hl.vds.filter_variants(vds_f, t, keep=True)
vds_f = hl.vds.filter_samples(vds_f, user_ids, keep = True, remove_dead_alleles = True)
mt = vds_f.variant_data.annotate_entries(
    AD = hl.vds.local_to_global(vds_f.variant_data.LAD, 
                                vds_f.variant_data.LA, 
                                n_alleles = hl.len(vds_f.variant_data.alleles), 
                                fill_value = 0, 
                                number = 'R'))
mt = mt.annotate_entries(GT = hl.vds.lgt_to_gt(mt.LGT, mt.LA))
mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds_f.reference_data, mt))

mt.GT.export("genotypes.tsv")
subprocess.run("hadoop fs -get genotypes.tsv", shell=True)
hl.stop()
