"""Extract methylation over insertion sites from a haplotagged PacBio modBAM.

For each insertion locus (BED), find primary reads that carry a large insertion
(an ``I`` CIGAR op) anchored at that locus, then report the inserted sequence,
the read's haplotype, and the 5mC/5hmC calls that fall inside the insertion --
both as a per-base string and as summary statistics. One output row per
(locus, read).
"""

import argparse
import csv

import pysam

# CIGAR ops that consume the query / reference (SAM spec 1.4.6).
# M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8
CONSUME_QUERY = {0, 1, 4, 7, 8}  # M, I, S, =, X
CONSUME_REF = {0, 2, 3, 7, 8}    # M, D, N, =, X

# Map single-letter and ChEBI integer mod codes to friendly labels.
MOD_CODE_ALIASES = {
    "m": "5mC", 27551: "5mC",
    "h": "5hmC", 76792: "5hmC",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True,
                        help="Indexed, haplotagged modBAM (.bai co-located).")
    parser.add_argument("--loci-bed", required=True,
                        help="BED of insertion loci: chrom, start, end, name, "
                             "min_len, max_len. min_len/max_len bound the "
                             "expected insertion length (bp); the accepted range "
                             "is widened by --len-tolerance on each side.")
    parser.add_argument("--output", default="insertion_methylation.tsv")
    parser.add_argument("--len-tolerance", type=int, default=100,
                        help="bp added on each side of each locus's "
                             "[min_len, max_len]; an insertion of length L is "
                             "accepted when min_len - tol <= L <= max_len + tol.")
    parser.add_argument("--anchor-pad", type=int, default=100,
                        help="bp window around the BED interval within which an "
                             "insertion anchor is accepted.")
    parser.add_argument("--mod-codes", default="m,h",
                        help="Comma-separated mod codes to keep (e.g. 'm,h').")
    parser.add_argument("--prob-threshold", type=float, default=0.5,
                        help="Probability >= threshold counts as methylated.")
    return parser.parse_args()


def read_bed(path):
    """Yield (chrom, start, end, name, min_len, max_len), skipping headers.

    Expects at least 6 columns: chrom, start, end, name, min_len, max_len.
    min_len/max_len bound the expected insertion length in bp.
    """
    with open(path) as inf:
        for i, line in enumerate(inf, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if len(fields) < 6:
                raise ValueError(
                    f"{path} line {i}: expected >=6 columns (chrom, start, end, "
                    f"name, min_len, max_len), got {len(fields)}: {line!r}")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if fields[3] else f"{chrom}:{start}-{end}"
            min_len, max_len = int(fields[4]), int(fields[5])
            yield chrom, start, end, name, min_len, max_len


def find_insertions(read, min_len, max_len):
    """Yield (q_start, ins_len, ref_anchor) for I ops in the length range.

    q_start indexes read.query_sequence; ref_anchor is the 0-based reference
    coordinate immediately before the insertion.
    """
    qpos = 0
    rpos = read.reference_start
    for op, length in read.cigartuples:
        if op == pysam.CINS and min_len <= length <= max_len:
            yield qpos, length, rpos - 1
        if op in CONSUME_QUERY:
            qpos += length
        if op in CONSUME_REF:
            rpos += length


def mods_in_span(read, q0, q1, keep_codes):
    """Return mod calls whose query index falls in [q0, q1).

    Uses read.modified_bases, whose positions index query_sequence -- the same
    frame as the CIGAR-derived offsets. Probability = (qual + 0.5) / 256;
    qual == -1 means no call and is dropped.
    """
    calls = []
    modified = read.modified_bases
    if not modified:
        return calls
    for (_canon, _strand, mod_code), positions in modified.items():
        if mod_code not in keep_codes:
            continue
        label = MOD_CODE_ALIASES.get(mod_code, str(mod_code))
        for qpos, qual in positions:
            if qual < 0 or not (q0 <= qpos < q1):
                continue
            calls.append((qpos - q0, label, (qual + 0.5) / 256.0))
    calls.sort()
    return calls


def summarize(calls, label, threshold):
    """Return (num_calls, num_methylated, mean_prob) for one mod label."""
    probs = [prob for _off, lab, prob in calls if lab == label]
    if not probs:
        return 0, 0, ""
    num_meth = sum(1 for p in probs if p >= threshold)
    mean_prob = round(sum(probs) / len(probs), 4)
    return len(probs), num_meth, mean_prob


COLUMNS = [
    "locus_name", "chrom", "insertion_ref_pos", "read_name", "haplotype",
    "phase_set", "strand", "insertion_length", "num_C",
    "num_5mC_calls", "num_5mC_methylated", "mean_5mC_prob",
    "num_5hmC_calls", "num_5hmC_methylated", "mean_5hmC_prob",
    "methylation_string", "insertion_sequence",
]


def main():
    args = parse_args()
    keep_codes = set()
    for code in args.mod_codes.split(","):
        code = code.strip()
        if not code:
            continue
        keep_codes.add(int(code) if code.isdigit() else code)

    bam = pysam.AlignmentFile(args.bam, "rb")

    with open(args.output, "w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(COLUMNS)

        for chrom, start, end, name, min_len, max_len in read_bed(args.loci_bed):
            lo = start - args.anchor_pad
            hi = end + args.anchor_pad
            min_ins_len = max(1, min_len - args.len_tolerance)
            max_ins_len = max_len + args.len_tolerance
            for read in bam.fetch(chrom, max(0, lo), hi):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                # Largest qualifying insertion anchored within the window.
                best = None
                for q_start, ins_len, anchor in find_insertions(
                        read, min_ins_len, max_ins_len):
                    if lo <= anchor <= hi and (best is None or ins_len > best[1]):
                        best = (q_start, ins_len, anchor)
                if best is None:
                    continue
                q_start, ins_len, anchor = best

                ins_seq = read.query_sequence[q_start:q_start + ins_len]
                haplotype = read.get_tag("HP") if read.has_tag("HP") else "NA"
                phase_set = read.get_tag("PS") if read.has_tag("PS") else "NA"
                strand = "-" if read.is_reverse else "+"

                calls = mods_in_span(read, q_start, q_start + ins_len, keep_codes)
                meth_string = ",".join(
                    f"{lab}:{off}:{prob:.3f}" for off, lab, prob in calls)

                n_5mc, m_5mc, mean_5mc = summarize(calls, "5mC", args.prob_threshold)
                n_5hmc, m_5hmc, mean_5hmc = summarize(calls, "5hmC", args.prob_threshold)

                writer.writerow([
                    name, chrom, anchor, read.query_name, haplotype, phase_set,
                    strand, ins_len, ins_seq.count("C"),
                    n_5mc, m_5mc, mean_5mc,
                    n_5hmc, m_5hmc, mean_5hmc,
                    meth_string, ins_seq,
                ])

    bam.close()


if __name__ == "__main__":
    main()
