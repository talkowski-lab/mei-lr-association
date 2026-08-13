"""Match assembly haplotypes to phased-VCF haplotypes over a region.

Assemblies are aligned to the reference (one BAM per haplotype), but the
haplotype numbering in the assembly has no relationship to the "0|1" phasing
in a VCF called independently from the same sample. This script figures out
the correspondence: for every phased, heterozygous, biallelic SNV the sample
has in the region, it reads the base each assembly haplotype carries at that
reference position (via pileup, no assembly variant calling needed) and
compares it to the two VCF-phased alleles. It reports, per phase set, the
percent concordance of each assembly haplotype against each VCF haplotype and
which pairing (hap1<->A/hap2<->B, or the swap) fits best.

VCF phase is only guaranteed to be consistent within a phase set (FORMAT/PS):
relative phase across two different PS blocks is arbitrary, so the region can
span a genuine phase switch. Concordance is therefore reported per phase set
(plus a pooled "ALL" row) rather than assuming one global phase.
"""

import argparse
import csv
import os
import subprocess

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hap1-bam", required=True,
                        help="Indexed BAM of haplotype-1 contigs aligned to the reference. "
                             "Local path, https://, or gs:// (with its .bai alongside it).")
    parser.add_argument("--hap2-bam", required=True,
                        help="Indexed BAM of haplotype-2 contigs aligned to the reference. "
                             "Local path, https://, or gs:// (with its .bai alongside it).")
    parser.add_argument("--vcf", required=True,
                        help="Bgzipped, tabix-indexed, phased multi-sample VCF. Local path, "
                             "https://, or gs:// (with its .tbi alongside it); for gs:// on a "
                             "private bucket, authenticate first (e.g. `gcloud auth login`) "
                             "or set GCS_OAUTH_TOKEN yourself.")
    parser.add_argument("--sample-id", required=True,
                        help="Sample column in the VCF to use as ground truth.")
    parser.add_argument("--region", required=True,
                        help="Region to compare, e.g. chr1:1000000-1010000 (1-based, inclusive).")
    parser.add_argument("--prefix", required=True,
                        help="Output basename prefix; writes '<prefix>.persite.tsv' and "
                             "'<prefix>.summary.tsv'.")
    parser.add_argument("--min-mapq", type=int, default=5,
                        help="Minimum mapping quality for an assembly alignment to be used.")
    parser.add_argument("--include-filtered", action="store_true",
                        help="Also use VCF records whose FILTER is set to something other "
                             "than PASS/missing (default: skip them).")
    return parser.parse_args()


def ensure_gcs_auth(paths):
    """Populate GCS_OAUTH_TOKEN from the local gcloud login, if needed.

    htslib reads gs:// paths itself (no code-level streaming needed -- it does
    ranged HTTP reads against the .tbi/.bai index, so it never downloads a
    whole cohort VCF just to read one region), but it only authenticates via
    the GCS_OAUTH_TOKEN env var. Without it, only public buckets are reachable.
    """
    if os.environ.get("GCS_OAUTH_TOKEN") or not any(p.startswith("gs://") for p in paths):
        return
    try:
        token = subprocess.run(
            ["gcloud", "auth", "print-access-token"],
            check=True, capture_output=True, text=True, timeout=30,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return
    if token:
        os.environ["GCS_OAUTH_TOKEN"] = token


def het_phased_snvs(vcf_path, sample_id, region, include_filtered):
    """Yield (chrom, pos0, ref, alt, hapA_allele, hapB_allele, phase_set).

    Restricted to biallelic SNVs where the sample's genotype is phased and
    heterozygous; hapA/hapB follow GT order (0|1 -> hapA=REF, hapB=ALT).
    """
    vcf = pysam.VariantFile(vcf_path)
    if sample_id not in vcf.header.samples:
        raise SystemExit(
            f"sample {sample_id!r} not found in {vcf_path} "
            f"(available: {', '.join(vcf.header.samples)})")
    # Restricts FORMAT parsing to this one sample instead of every column in
    # the cohort VCF -- the pysam/htslib equivalent of `bcftools view -s`.
    vcf.subset_samples([sample_id])

    for rec in vcf.fetch(region=region):
        if not include_filtered:
            filters = list(rec.filter.keys())
            if filters and filters != ["PASS"]:
                continue
        if len(rec.alleles) != 2:
            continue
        ref, alt = rec.alleles
        if len(ref) != 1 or len(alt) != 1:
            continue

        sample = rec.samples[sample_id]
        gt = sample["GT"]
        if gt[0] is None or gt[1] is None or not sample.phased:
            continue
        if gt[0] == gt[1]:
            continue

        alleles = (ref, alt)
        hapA_allele = alleles[gt[0]]
        hapB_allele = alleles[gt[1]]
        phase_set = sample["PS"] if "PS" in sample else "NA"
        yield rec.chrom, rec.pos - 1, ref, alt, hapA_allele, hapB_allele, str(phase_set)


def base_at(bam, chrom, pos0, min_mapq):
    """Return the base an assembly BAM carries at pos0, "AMBIG", or None."""
    for column in bam.pileup(
            chrom, pos0, pos0 + 1, truncate=True,
            min_base_quality=0, min_mapping_quality=min_mapq):
        if column.reference_pos != pos0:
            continue
        bases = set()
        for pileupread in column.pileups:
            aln = pileupread.alignment
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if pileupread.is_del or pileupread.is_refskip:
                continue
            if pileupread.query_position is None:
                continue
            bases.add(aln.query_sequence[pileupread.query_position].upper())
        if not bases:
            return None
        if len(bases) > 1:
            return "AMBIG"
        return bases.pop()
    return None


PERSITE_COLUMNS = [
    "chrom", "pos", "ref", "alt", "phase_set",
    "vcf_hapA_allele", "vcf_hapB_allele", "asm_hap1_base", "asm_hap2_base",
]

PAIRS = [("hap1", "A"), ("hap1", "B"), ("hap2", "A"), ("hap2", "B")]


def new_tally():
    return {pair: [0, 0] for pair in PAIRS}  # pair -> [matches, informative]


def update_tally(tally, asm_bases, vcf_alleles):
    for pair in PAIRS:
        asm_hap, vcf_hap = pair
        asm_base = asm_bases[asm_hap]
        if asm_base is None or asm_base == "AMBIG":
            continue
        tally[pair][1] += 1
        if asm_base == vcf_alleles[vcf_hap]:
            tally[pair][0] += 1


def pct(matches, informative):
    return round(100.0 * matches / informative, 2) if informative else None


def best_assignment(tally):
    n1a, d1a = tally[("hap1", "A")]
    n2b, d2b = tally[("hap2", "B")]
    n1b, d1b = tally[("hap1", "B")]
    n2a, d2a = tally[("hap2", "A")]

    straight = (n1a + n2b) / (d1a + d2b) if (d1a + d2b) else None
    swapped = (n1b + n2a) / (d1b + d2a) if (d1b + d2a) else None

    if straight is None and swapped is None:
        return "NA", None
    if swapped is None or (straight is not None and straight >= swapped):
        return "hap1=A,hap2=B", round(straight * 100, 2)
    return "hap1=B,hap2=A", round(swapped * 100, 2)


def write_summary(path, tallies):
    with open(path, "w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow([
            "phase_set", "n_hap1_vs_A", "hap1_vs_A_pct",
            "n_hap1_vs_B", "hap1_vs_B_pct",
            "n_hap2_vs_A", "hap2_vs_A_pct",
            "n_hap2_vs_B", "hap2_vs_B_pct",
            "best_assignment", "best_assignment_pct",
        ])
        for phase_set in sorted(tallies, key=lambda k: (k != "ALL", k)):
            tally = tallies[phase_set]
            assignment, assignment_pct = best_assignment(tally)
            row = [phase_set]
            for pair in PAIRS:
                matches, informative = tally[pair]
                row += [informative, pct(matches, informative)]
            row += [assignment, assignment_pct]
            writer.writerow(row)


def main():
    args = parse_args()
    ensure_gcs_auth([args.hap1_bam, args.hap2_bam, args.vcf])

    hap1_bam = pysam.AlignmentFile(args.hap1_bam, "rb")
    hap2_bam = pysam.AlignmentFile(args.hap2_bam, "rb")

    tallies = {"ALL": new_tally()}

    persite_path = f"{args.prefix}.persite.tsv"
    with open(persite_path, "w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t")
        writer.writerow(PERSITE_COLUMNS)

        for chrom, pos0, ref, alt, hapA, hapB, phase_set in het_phased_snvs(
                args.vcf, args.sample_id, args.region, args.include_filtered):
            hap1_base = base_at(hap1_bam, chrom, pos0, args.min_mapq)
            hap2_base = base_at(hap2_bam, chrom, pos0, args.min_mapq)

            writer.writerow([
                chrom, pos0 + 1, ref, alt, phase_set,
                hapA, hapB, hap1_base or "NA", hap2_base or "NA",
            ])

            asm_bases = {"hap1": hap1_base, "hap2": hap2_base}
            vcf_alleles = {"A": hapA, "B": hapB}
            update_tally(tallies["ALL"], asm_bases, vcf_alleles)
            tallies.setdefault(phase_set, new_tally())
            update_tally(tallies[phase_set], asm_bases, vcf_alleles)

    hap1_bam.close()
    hap2_bam.close()

    write_summary(f"{args.prefix}.summary.tsv", tallies)


if __name__ == "__main__":
    main()
