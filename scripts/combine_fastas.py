import sys
from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(description="Combine multiple FASTA files into one.")
parser.add_argument("input_files", nargs="+", help="Input FASTA files to combine.")
args = parser.parse_args()


records = []
for f in args.input_files:
    with open(f, "r") as handle:
        records.extend(list(SeqIO.parse(handle, "fasta")))

records = sorted(records, key=lambda x: x.id)
# Get rid of short sequences
records = [r for r in records if len(r.seq) >= 100]
SeqIO.write(records, sys.stdout, "fasta")
