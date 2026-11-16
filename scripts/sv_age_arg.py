#!/usr/bin/env python3
"""Date a structural variant from a local ARG.

Infers a tree sequence (tsinfer) from the phased SNVs flanking an SV, with the
SV's own span excluded, adds the singleton SNVs from the unphased callset as
non-inference sites, dates it (tsdate), then maps the SV onto the trees by
parsimony and reports its age, the mutations needed to explain it, and the
individuals whose calls disagree with the genealogy.
"""

import argparse
import gzip
import json
import logging
import os
import shutil
import sys
import time

import cyvcf2
import numpy as np
import pandas as pd
import pysam
import tsdate
import tsinfer
import tskit
import tszip

log = logging.getLogger("sv_age_arg")
NUCS = frozenset("ACGT")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    io = p.add_argument_group("inputs")
    io.add_argument("--phased-vcf", required=True,
                    help="Phased MAC>=2 VCF covering the flank region (biallelic SNVs used).")
    io.add_argument("--singleton-vcf",
                    help="Full-callset VCF covering the flank region; SNVs whose derived allele "
                         "is carried once among the phased VCF's samples are added for dating.")
    io.add_argument("--sv-phased-vcf", help="Records near the SV from the phased VCF.")
    io.add_argument("--sv-full-vcf", help="Records near the SV from the full VCF.")
    io.add_argument("--sv-ids", required=True,
                    help="Comma-separated SV entry IDs: the VCF ID, CHROM:POS:REF:ALT or "
                         "CHROM:POS_REF_ALT. A haplotype carries the SV if it carries any entry.")
    io.add_argument("--chrom", required=True)
    io.add_argument("--sv-start", type=int, required=True, help="1-based SV start.")
    io.add_argument("--sv-end", type=int, required=True, help="1-based SV end (inclusive).")
    io.add_argument("--flank", type=int, default=500_000)
    io.add_argument("--exclude-pad", type=int, default=500,
                    help="Also drop SNVs this close to the SV span.")
    io.add_argument("--prefix", required=True)

    anc = p.add_argument_group("polarisation")
    anc.add_argument("--ancestral-fasta",
                     help="Ancestral-allele FASTA (e.g. Ensembl homo_sapiens_ancestor; "
                          "plain or bgzipped, whole-genome or one chromosome).")
    anc.add_argument("--ancestral-fasta-index")
    anc.add_argument("--ancestral-offset", type=int,
                     help="The FASTA holds one region, not whole chromosomes: its first "
                          "base is at this 1-based position on --chrom.")
    anc.add_argument("--ancestral-info-field",
                     help="INFO field holding the ancestral allele, used if no FASTA.")
    anc.add_argument("--ancestral-high-confidence-only", action="store_true",
                     help="Ignore lower-case (low-confidence) FASTA calls.")
    anc.add_argument("--sv-type", type=str.upper, choices=["INS", "DEL"], default="INS",
                     help="INS: REF is ancestral (the insertion is dated). DEL: ALT is "
                          "ancestral (the reference sequence is the derived insertion, e.g. a "
                          "reference SVA missing in some people).")
    anc.add_argument("--sv-ancestral", choices=["ref", "alt", "auto"],
                     help="Override --sv-type's polarity; auto = whichever parsimony prefers.")

    model = p.add_argument_group("model")
    model.add_argument("--genetic-map",
                       help="Genetic map: chrom/pos/rate/cM (HapMap/SHAPEIT style), PLINK .map, "
                            "or pos/rate/cM for one chromosome. Header optional, .gz ok.")
    model.add_argument("--recombination-rate", type=float, default=1e-8,
                       help="Per-bp per-generation rate used when no map is given.")
    model.add_argument("--mismatch-ratio", type=float, default=1.0)
    model.add_argument("--mutation-rate", type=float, default=1.29e-8)
    model.add_argument("--generation-time", type=float, default=28.0)
    model.add_argument("--singletons-phased", action="store_true",
                       help="Treat singleton phase as known (only if the singleton VCF is phased).")
    model.add_argument("--threads", type=int, default=4)

    qc = p.add_argument_group("qc")
    qc.add_argument("--scan-bp", type=int, default=20_000,
                    help="Also map the SV on every tree this close to it.")
    qc.add_argument("--background-bp", type=int, default=100_000,
                    help="Window for the SNV parsimony-cost background.")
    qc.add_argument("--background-max-sites", type=int, default=200)
    qc.add_argument("--max-error-clade-haps", type=int, default=2,
                    help="An extra mutation explaining at most this many haplotypes is treated "
                         "as a genotyping/phasing error.")
    qc.add_argument("--recurrent-min-haps", type=int, default=5,
                    help="An extra forward mutation explaining at least this many carrier "
                         "haplotypes is treated as a possible recurrent event.")
    qc.add_argument("--min-flank-sites", type=int, default=10,
                    help="Flag low site density when fewer inference sites than this lie within "
                         "10 kb on either side of the SV.")
    qc.add_argument("--max-phase-tests", type=int, default=2000)
    qc.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def norm_chrom(c):
    c = str(c)
    if c.lower().startswith("chr"):
        c = c[3:]
    c = c.upper()
    return {"23": "X", "24": "Y", "M": "MT"}.get(c, c)


# ---------------------------------------------------------------- ancestral alleles

def open_fasta(path, index=None):
    local = os.path.basename(path)
    if os.path.abspath(path) != os.path.abspath(local) and not os.path.exists(local):
        os.symlink(os.path.abspath(path), local)
    for ext, src in ((".fai", index), (".fai", path + ".fai"), (".gzi", path + ".gzi")):
        if src and os.path.exists(src) and not os.path.exists(local + ext):
            os.symlink(os.path.abspath(src), local + ext)
    try:
        return pysam.FastaFile(local)
    except (OSError, ValueError):
        if not path.endswith(".gz"):
            raise
        plain = local[:-3]
        log.info("FASTA is gzip but not bgzip; decompressing to %s", plain)
        with gzip.open(path, "rt") as fi, open(plain, "w") as fo:
            shutil.copyfileobj(fi, fo)
        return pysam.FastaFile(plain)


def find_contig(fa, chrom):
    want = norm_chrom(chrom)
    for ref in fa.references:
        if ref == chrom or norm_chrom(ref) == want:
            return ref
    # Ensembl ancestral headers: ANCESTOR_for_chromosome:GRCh38:1:1:248956422:1
    for ref in fa.references:
        parts = ref.split(":")
        if len(parts) >= 3 and parts[0].startswith("ANCESTOR") and norm_chrom(parts[2]) == want:
            return ref
    raise SystemExit(f"Contig {chrom} not found in ancestral FASTA "
                     f"(first contigs: {', '.join(fa.references[:5])})")


class Ancestral:
    def __init__(self, args, start, end):
        self.field = None
        self.seq = None
        self.high_conf_only = args.ancestral_high_confidence_only
        if args.ancestral_fasta:
            fa = open_fasta(args.ancestral_fasta, args.ancestral_fasta_index)
            if args.ancestral_offset:
                contig = fa.references[0]
                self.offset = args.ancestral_offset
                self.seq = fa.fetch(contig)
            else:
                contig = find_contig(fa, args.chrom)
                self.offset = start
                self.seq = fa.fetch(contig, start - 1, end)
            self.source = f"fasta:{os.path.basename(args.ancestral_fasta)}:{contig}"
        elif args.ancestral_info_field:
            self.field = args.ancestral_info_field
            self.source = f"info:{self.field}"
        else:
            self.source = "reference"
            log.warning("No ancestral alleles given: REF is treated as ancestral, which "
                        "mis-polarises a few percent of sites and biases dates.")

    def get(self, v):
        if self.seq is not None:
            i = v.POS - self.offset
            if i < 0 or i >= len(self.seq):
                return None
            a = self.seq[i]
            if self.high_conf_only and a.islower():
                return None
            a = a.upper()
        elif self.field:
            a = v.INFO.get(self.field)
            if a is None:
                return None
            a = str(a).split(",")[0].strip().upper()
        else:
            a = v.REF.upper()
        return a if a in NUCS else None


# ---------------------------------------------------------------- genetic map

def read_genetic_map(path, chrom, seq_len, fallback_rate):
    opener = gzip.open if path.endswith(".gz") else open
    want = norm_chrom(chrom)
    pos, cm = [], []
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if not f or f[0].startswith("#"):
                continue
            try:
                if len(f) >= 4 and not _isnum(f[1]):          # PLINK: chrom id cM bp
                    c, p_, m = f[0], int(float(f[3])), float(f[2])
                elif len(f) >= 4:                              # chrom pos rate cM
                    c, p_, m = f[0], int(float(f[1])), float(f[3])
                elif len(f) == 3:                              # pos rate cM
                    c, p_, m = chrom, int(float(f[0])), float(f[2])
                else:
                    continue
            except ValueError:
                continue                                       # header
            if norm_chrom(c) == want:
                pos.append(p_)
                cm.append(m)
    if len(pos) < 2:
        raise SystemExit(f"Genetic map {path} has fewer than two points for {chrom}")
    order = np.argsort(pos, kind="stable")
    pos, cm = np.asarray(pos)[order], np.asarray(cm)[order]
    keep = np.concatenate([[True], np.diff(pos) > 0])
    pos, cm = pos[keep], cm[keep]
    seg_rate = np.diff(cm) / 100 / np.diff(pos)
    seg_rate = np.maximum(seg_rate, 0)
    inside = (pos > 0) & (pos < seq_len)
    mean_rate = (cm[-1] - cm[0]) / 100 / (pos[-1] - pos[0]) if pos[-1] > pos[0] else fallback_rate
    breaks = [0.0]
    rates = []
    # rate for [0, first inside point): the segment ending there, else the mean
    idx = np.flatnonzero(inside)
    if len(idx) == 0:
        log.warning("No genetic-map points inside the region; using the chromosome mean rate")
        return tskit.RateMap(position=[0, seq_len], rate=[max(mean_rate, 1e-12)]), "map-mean"
    first = idx[0]
    rates.append(seg_rate[first - 1] if first > 0 else mean_rate)
    for j in idx:
        breaks.append(float(pos[j]))
        rates.append(seg_rate[j] if j < len(seg_rate) else mean_rate)
    breaks.append(float(seq_len))
    rates = np.maximum(np.asarray(rates, float), 1e-12)
    return tskit.RateMap(position=breaks, rate=rates), f"map:{os.path.basename(path)}"


def _isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------- VCF parsing

def gt_matrix(v):
    """(n_samples, 2) allele matrix (-1 missing, -2 absent second copy) and phased flags."""
    a = v.genotype.array()
    g = a[:, :-1]
    if g.shape[1] == 1:
        g = np.column_stack([g, np.full(len(g), -2, dtype=g.dtype)])
    elif g.shape[1] > 2:
        raise SystemExit(f"Ploidy > 2 at {v.CHROM}:{v.POS}; not supported")
    return g, a[:, -1] == 1


def is_biallelic_snv(v):
    return (len(v.ALT) == 1 and len(v.REF) == 1 and v.REF.upper() in NUCS
            and v.ALT[0].upper() in NUCS)


def sample_layout(path, n_probe=500):
    vcf = cyvcf2.VCF(path)
    samples = list(vcf.samples)
    ploidy = np.zeros(len(samples), dtype=int)
    for k, v in enumerate(vcf):
        g, _ = gt_matrix(v)
        ploidy = np.maximum(ploidy, (g != -2).sum(1))
        if k >= n_probe:
            break
    vcf.close()
    ploidy[ploidy == 0] = 2
    hap_cols = np.array([i * 2 + j for i in range(len(samples)) for j in range(ploidy[i])])
    hap_owner = np.repeat(np.arange(len(samples)), ploidy)
    return samples, ploidy, hap_cols, hap_owner


class Region:
    def __init__(self, args):
        self.chrom = args.chrom
        self.sv_start, self.sv_end = sorted((args.sv_start, args.sv_end))
        self.start = max(1, self.sv_start - args.flank)
        self.end = self.sv_end + args.flank
        self.excl_lo = self.sv_start - args.exclude_pad
        self.excl_hi = self.sv_end + args.exclude_pad
        self.seq_len = self.end + 1

    def keep(self, pos):
        return self.start <= pos <= self.end and not (self.excl_lo <= pos <= self.excl_hi)


def read_singletons(path, region, anc, samples, ploidy, stats):
    """Sites whose derived allele occurs exactly once among the phased VCF's samples."""
    out = {}
    if not path:
        return out
    vcf = cyvcf2.VCF(path)
    index = {s: i for i, s in enumerate(samples)}
    f_idx = np.array([index.get(s, -1) for s in vcf.samples])
    shared_rows = np.flatnonzero(f_idx >= 0)
    stats["singleton_vcf_samples_not_in_phased"] = int((f_idx < 0).sum())
    stats["phased_samples_not_in_singleton_vcf"] = len(samples) - len(shared_rows)
    hap_start = np.concatenate([[0], np.cumsum(ploidy)[:-1]])
    absent = np.ones(int(ploidy.sum()), dtype=bool)
    for r in shared_rows:
        p = f_idx[r]
        absent[hap_start[p]:hap_start[p] + ploidy[p]] = False
    absent_haps = np.flatnonzero(absent)
    for v in vcf:
        if not is_biallelic_snv(v) or not region.keep(v.POS):
            continue
        a = anc.get(v)
        ref, alt = v.REF.upper(), v.ALT[0].upper()
        if a not in (ref, alt):
            continue
        g, _ = gt_matrix(v)
        g = g[shared_rows]
        if a == alt:
            g = np.where(g >= 0, 1 - g, g)
        if int((g == 1).sum()) != 1:
            continue
        row, col = map(int, np.argwhere(g == 1)[0])
        p = f_idx[shared_rows[row]]
        carrier = hap_start[p] + min(col, ploidy[p] - 1)
        miss = []
        for rr in np.flatnonzero((g == -1).any(1)):
            pp = f_idx[shared_rows[rr]]
            for j in range(ploidy[pp]):
                if g[rr, j] == -1:
                    miss.append(hap_start[pp] + j)
        der = alt if a == ref else ref
        out[v.POS] = (a, der, carrier, np.asarray(miss, dtype=np.int64))
    vcf.close()
    out["__absent__"] = absent_haps
    return out


def build_sample_data(args, region, anc, samples, ploidy, hap_cols, singletons, stats):
    n_haps = int(ploidy.sum())
    sd = tsinfer.SampleData(path=f"{args.prefix}.samples", sequence_length=region.seq_len,
                            num_flush_threads=2)
    for s, k in zip(samples, ploidy):
        sd.add_individual(ploidy=int(k), metadata={"sample_id": s})

    absent = singletons.pop("__absent__", np.array([], dtype=np.int64))
    single_pos = sorted(singletons)
    si = 0
    counters = dict.fromkeys(
        ["records", "not_biallelic_snv", "outside_region_or_in_sv", "no_ancestral",
         "ancestral_not_ref_or_alt", "monomorphic", "duplicate_position", "inference_sites",
         "singleton_sites", "unphased_het_calls_masked", "ancestral_is_alt"], 0)
    positions = []
    derived = []

    def add_singleton(pos):
        a, der, carrier, miss = singletons[pos]
        h = np.zeros(n_haps, dtype=np.int8)
        h[miss] = -1
        h[absent] = -1
        h[carrier] = 1
        sd.add_site(pos, h, alleles=[a, der])
        counters["singleton_sites"] += 1

    pending = None
    dup_pos = -1

    def flush(pos_limit):
        nonlocal si
        while si < len(single_pos) and single_pos[si] < pos_limit:
            add_singleton(single_pos[si])
            si += 1

    def commit(site):
        nonlocal si
        pos, h, alleles = site
        flush(pos)
        if si < len(single_pos) and single_pos[si] == pos:
            si += 1
        sd.add_site(pos, h, alleles=alleles)
        positions.append(pos)
        derived.append(int((h == 1).sum()))
        counters["inference_sites"] += 1

    vcf = cyvcf2.VCF(args.phased_vcf)
    if list(vcf.samples) != samples:
        raise SystemExit("Phased VCF sample order changed between reads")
    for v in vcf:
        counters["records"] += 1
        if not is_biallelic_snv(v):
            counters["not_biallelic_snv"] += 1
            continue
        if not region.keep(v.POS):
            counters["outside_region_or_in_sv"] += 1
            continue
        if pending is not None and v.POS == pending[0]:
            counters["duplicate_position"] += 2
            dup_pos = v.POS
            pending = None
            continue
        if v.POS == dup_pos:
            counters["duplicate_position"] += 1
            continue
        if pending is not None:
            commit(pending)
            pending = None
        a = anc.get(v)
        ref, alt = v.REF.upper(), v.ALT[0].upper()
        if a is None:
            counters["no_ancestral"] += 1
            continue
        if a not in (ref, alt):
            counters["ancestral_not_ref_or_alt"] += 1
            continue
        g, phased = gt_matrix(v)
        het = (g[:, 0] >= 0) & (g[:, 1] >= 0) & (g[:, 0] != g[:, 1]) & ~phased
        if het.any():
            counters["unphased_het_calls_masked"] += int(het.sum())
            g = g.copy()
            g[het] = -1
        h = g.reshape(-1)[hap_cols].astype(np.int8)
        h[h < -1] = -1
        if a == alt:
            counters["ancestral_is_alt"] += 1
            h = np.where(h >= 0, 1 - h, h).astype(np.int8)
        called = h >= 0
        dc = int((h == 1).sum())
        if dc == 0 or dc == int(called.sum()):
            counters["monomorphic"] += 1
            continue
        pending = (v.POS, h, [a, alt if a == ref else ref])
    vcf.close()
    if pending is not None:
        commit(pending)
    flush(region.seq_len + 1)
    sd.finalise()
    stats.update({f"phased_{k}" if k not in ("inference_sites", "singleton_sites") else k: val
                  for k, val in counters.items()})
    return sd, np.asarray(positions), np.asarray(derived)


# ---------------------------------------------------------------- SV genotypes

def id_keys(v):
    keys = set()
    if v.ID:
        keys.update(x for x in v.ID.split(";") if x and x != ".")
    for alt in v.ALT:
        keys.add(f"{v.CHROM}:{v.POS}:{v.REF}:{alt}")
        keys.add(f"{v.CHROM}:{v.POS}_{v.REF}_{alt}")
    return keys


def collect_sv_entries(path, wanted):
    found = {}
    nearby = []
    if not path or not os.path.exists(path):
        return found, nearby
    vcf = cyvcf2.VCF(path)
    for v in vcf:
        keys = id_keys(v)
        nearby.append(v.ID or f"{v.CHROM}:{v.POS}")
        for w in wanted:
            if w in keys and w not in found:
                g, phased = gt_matrix(v)
                found[w] = (list(vcf.samples), g.copy(), phased.copy(), v.POS)
    vcf.close()
    return found, nearby


def sv_haplotypes(args, samples, ploidy):
    """Per-haplotype SV state (1 carrier, 0 not, -1 missing) and individuals whose
    carrier haplotype is unknown because the call was unphased."""
    wanted = [s.strip() for s in args.sv_ids.split(",") if s.strip()]
    ph, near_p = collect_sv_entries(args.sv_phased_vcf, wanted)
    fu, near_f = collect_sv_entries(args.sv_full_vcf, [w for w in wanted if w not in ph])
    missing = [w for w in wanted if w not in ph and w not in fu]
    if len(missing) == len(wanted):
        near = sorted(set(near_p + near_f))[:30]
        raise SystemExit(f"None of the SV entries {wanted} were found near "
                         f"{args.chrom}:{args.sv_start}-{args.sv_end}. Records there: {near}")
    if missing:
        log.warning("SV entries not found and ignored: %s", missing)

    n = len(samples)
    index = {s: i for i, s in enumerate(samples)}
    hap_start = np.concatenate([[0], np.cumsum(ploidy)[:-1]])
    carrier = np.zeros((n, 2), dtype=bool)
    called = np.zeros((n, 2), dtype=bool)
    unphased = np.zeros(n, dtype=bool)
    entry_rows = []
    for source, entries in (("phased", ph), ("full", fu)):
        for w, (vs, g, phased, pos) in entries.items():
            rows = np.array([index.get(s, -1) for s in vs])
            ok = rows >= 0
            r = rows[ok]
            g, phased = g[ok], phased[ok]
            if source == "full":
                # the full callset's phase (e.g. per-sample read-backed) does not share the
                # phased VCF's haplotype labels, so its "|" is ignored
                phased = np.zeros_like(phased)
            car = g > 0
            called[r] |= g >= 0
            het = car.any(1) & ~car.all(1) & (g >= 0).all(1)
            unph = het & ~phased & ~(ploidy[r] == 1)
            # an unphased het never picks a haplotype itself: an existing (phased) carrier
            # keeps its placement, otherwise it starts on haplotype 0 and the tree decides
            known = carrier[r].any(1)
            car[unph] = False
            car[unph & ~known, 0] = True
            carrier[r] |= car
            unphased[r[unph & ~known]] = True
            entry_rows.append({"entry": w, "source": source, "pos": pos,
                               "carrier_individuals": int(car.any(1).sum()),
                               "unphased_het_carriers": int(unph.sum()),
                               "missing_individuals": int((g == -1).all(1).sum())})

    # a non-carrier is missing only when no entry calls it: a split MEI often has a sparse
    # record that is "./." for most samples, which must not erase their calls elsewhere
    state = np.zeros((n, 2), dtype=np.int8)
    state[carrier] = 1
    state[~carrier & ~called] = -1
    # a carrier with an unphased call goes on its first haplotype until the tree places it
    for i in np.flatnonzero(unphased):
        if carrier[i].sum() == 1:
            state[i] = (1, 0)
        else:
            unphased[i] = False
    h = np.concatenate([state[i, :ploidy[i]] for i in range(n)])
    unphased_inds = np.flatnonzero(unphased)
    info = {"sv_ids": wanted, "found_phased": sorted(ph), "found_full": sorted(fu),
            "not_found": missing, "entries": entry_rows,
            "n_individuals_no_sv_call": int((~called.any(1)).sum())}
    return h, unphased_inds, hap_start, info


# ---------------------------------------------------------------- parsimony on a tree

def fitch(tree, geno, mode):
    alleles = ["0", "1"]
    if mode == "auto":
        a, trans = tree.map_mutations(geno, alleles)
    else:
        a, trans = tree.map_mutations(geno, alleles, ancestral_state="0" if mode == "ref" else "1")
    return int(a), trans


def cost(tree, geno, mode):
    return len(fitch(tree, geno, mode)[1])


def describe(tree, geno, anc, trans):
    """Assign each sample to the deepest mutation above it and pick the main mutation."""
    der = 1 - anc
    n = len(geno)
    depth = []
    for m in trans:
        d, p = 0, m.parent
        while p != tskit.NULL:
            d, p = d + 1, trans[p].parent
        depth.append(d)
    owner = np.full(n, -1, dtype=np.int64)
    under = []
    for k in np.argsort(depth, kind="stable"):
        s = np.fromiter(tree.samples(trans[k].node), dtype=np.int64)
        owner[s] = k
    for m in trans:
        under.append(np.fromiter(tree.samples(m.node), dtype=np.int64))
    rows = []
    main = None
    best = -1
    for k, m in enumerate(trans):
        mine = owner == k
        r = {"k": k, "node": m.node, "state": int(m.derived_state), "parent_mut": m.parent,
             "haps_total": len(under[k]), "haps_exclusive": int(mine.sum()),
             "derived_called": int((mine & (geno == der)).sum()),
             "ancestral_called": int((mine & (geno == anc)).sum()),
             "missing": int((mine & (geno < 0)).sum())}
        rows.append(r)
        if r["state"] == der and r["derived_called"] > best:
            best, main = r["derived_called"], k
    return owner, under, rows, main


def ancestors_of(trans, k):
    out = []
    p = trans[k].parent
    while p != tskit.NULL:
        out.append(p)
        p = trans[p].parent
    return out


def classify(rows, trans, main, anc, max_err, rec_min):
    der = 1 - anc
    for r in rows:
        k = r["k"]
        if k == main:
            r["role"] = "main"
            continue
        inside = main in ancestors_of(trans, k)
        forward = r["state"] == der
        r["inside_main"] = inside
        n_bad = r["derived_called"] if forward else r["ancestral_called"]
        if forward:
            if n_bad <= max_err:
                r["role"] = "isolated_carriers"
            elif n_bad >= rec_min and not inside:
                r["role"] = "recurrent_candidate"
            else:
                r["role"] = "ambiguous_carrier_clade"
        else:
            if n_bad <= max_err:
                r["role"] = "isolated_noncarriers"
            elif n_bad >= rec_min:
                r["role"] = "large_reversion"
            else:
                r["role"] = "ambiguous_noncarrier_clade"
    return rows


# ---------------------------------------------------------------- main analysis

def node_age(ts, node):
    t = float(ts.nodes_time[node])
    md = ts.node(node).metadata
    mn = vr = np.nan
    if isinstance(md, dict):
        mn, vr = md.get("mn", np.nan), md.get("vr", np.nan)
    return t, mn, vr


def branch_age(ts, tree, node):
    lo = float(ts.nodes_time[node])
    parent = tree.parent(node)
    hi = float(ts.nodes_time[parent]) if parent != tskit.NULL else np.nan
    mid = (lo + hi) / 2 if parent != tskit.NULL else np.nan
    return lo, hi, mid


def place_unphased(tree, geno, unphased_inds, hap_start, mode):
    """Put each unphased carrier's SV on whichever haplotype gives the lower cost."""
    placed = []
    base = cost(tree, geno, mode)
    for i in unphased_inds:
        a, b = hap_start[i], hap_start[i] + 1
        alt = geno.copy()
        alt[a], alt[b] = geno[b], geno[a]
        c = cost(tree, alt, mode)
        moved = c < base
        if moved:
            geno, base = alt, c
        placed.append((i, int(b if moved else a)))
    return geno, placed


def pick_holdout(inf_pos, inf_dc, sv_pos, sv_dc, args):
    """SNVs near the SV with a similar derived count, withheld from inference.

    Inference sites are fitted by construction, so they understate how often a
    single real mutation maps imperfectly; held-out sites are placed onto the
    trees afterwards by the same parsimony as the SV, which makes them a fair
    calibration."""
    if sv_dc < 2:
        return np.array([], dtype=float), (np.nan, np.nan)
    lo_c, hi_c = max(2, sv_dc / 2), sv_dc * 2
    d = np.abs(inf_pos - sv_pos)
    near = d <= args.background_bp
    # the trees the SV itself is mapped on keep every site
    cand = np.flatnonzero(near & (d > args.scan_bp) & (inf_dc >= lo_c) & (inf_dc <= hi_c))
    cap = min(args.background_max_sites, int(0.05 * near.sum()))
    if len(cand) > cap:
        cand = np.random.default_rng(args.seed).choice(cand, cap, replace=False)
    return np.sort(inf_pos[cand]).astype(float), (lo_c, hi_c)


def split_holdout(sd, holdout, path):
    """Inference SampleData without the held-out sites, plus their genotypes."""
    pos = sd.sites_position[:]
    is_out = np.isin(pos, holdout)
    held = [(float(pos[i]), np.asarray(sd.sites_genotypes[i], dtype=np.int32))
            for i in np.flatnonzero(is_out)]
    if not held:
        return sd, held
    sub = sd.subset(sites=np.flatnonzero(~is_out), path=path)
    return sub, held


def background_costs(ts, held, rec_min):
    """Parsimony cost of each held-out SNV (ancestral = 0), and whether it splits into
    >= 2 large forward clades, which is what a recurrent SV would look like."""
    costs, splits = [], []
    tree = ts.first()
    lo, hi = ts.first().interval.left, ts.last().interval.right
    for pos, g in held:
        if not (lo <= pos < hi):
            continue
        tree.seek(pos)
        _, tr = tree.map_mutations(g, ["0", "1"], ancestral_state=0)
        big = 0
        for m in tr:
            if m.derived_state == "1":
                s = np.fromiter(tree.samples(m.node), dtype=np.int64)
                big += int((g[s] == 1).sum()) >= rec_min
        costs.append(len(tr))
        splits.append(big >= 2)
    return np.asarray(costs, dtype=float), np.asarray(splits, dtype=bool)


def main(argv=None):
    args = parse_args(argv)
    if args.sv_ancestral is None:
        args.sv_ancestral = "ref" if args.sv_type == "INS" else "alt"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr),
                                  logging.FileHandler(f"{args.prefix}.log", mode="w")])
    t0 = time.time()
    region = Region(args)
    stats = {"chrom": args.chrom, "region_start": region.start, "region_end": region.end,
             "excluded_start": region.excl_lo, "excluded_end": region.excl_hi,
             "tsinfer_version": tsinfer.__version__, "tsdate_version": tsdate.__version__,
             "tskit_version": tskit.__version__}

    samples, ploidy, hap_cols, hap_owner = sample_layout(args.phased_vcf)
    n_haps = int(ploidy.sum())
    stats.update(n_individuals=len(samples), n_haplotypes=n_haps,
                 n_haploid_individuals=int((ploidy == 1).sum()))
    log.info("%d individuals, %d haplotypes", len(samples), n_haps)

    anc = Ancestral(args, region.start, region.end)
    stats["ancestral_source"] = anc.source

    singletons = read_singletons(args.singleton_vcf, region, anc, samples, ploidy, stats)
    log.info("%d singleton sites from the full callset", len(singletons) - ("__absent__" in singletons))

    geno, unphased_inds, hap_start, sv_info = sv_haplotypes(args, samples, ploidy)
    if (geno == 1).sum() == 0 or (geno == 0).sum() == 0:
        raise SystemExit("The SV is monomorphic among the phased samples")
    n_car, n_non = int((geno == 1).sum()), int((geno == 0).sum())
    sv_dc = {"ref": n_car, "alt": n_non}.get(args.sv_ancestral, min(n_car, n_non))

    sd, inf_pos, inf_dc = build_sample_data(args, region, anc, samples, ploidy, hap_cols,
                                            singletons, stats)
    holdout, (bg_lo, bg_hi) = pick_holdout(inf_pos, inf_dc, region.sv_start, sv_dc, args)
    sd, held = split_holdout(sd, holdout, f"{args.prefix}.infer.samples")
    stats["holdout_sites"] = len(held)
    inf_pos = np.setdiff1d(inf_pos, holdout.astype(inf_pos.dtype))
    if len(inf_pos) < 10:
        raise SystemExit(f"Only {len(inf_pos)} usable SNVs in the region; cannot build an ARG")
    stats["inference_sites_per_kb"] = round(len(inf_pos) / ((region.end - region.start) / 1e3), 3)
    left = inf_pos[inf_pos < region.sv_start]
    right = inf_pos[inf_pos > region.sv_end]
    stats["nearest_site_left_bp"] = int(region.sv_start - left.max()) if len(left) else -1
    stats["nearest_site_right_bp"] = int(right.min() - region.sv_end) if len(right) else -1
    stats["sites_within_10kb_left"] = int((left >= region.sv_start - 10_000).sum())
    stats["sites_within_10kb_right"] = int((right <= region.sv_end + 10_000).sum())
    if not len(left) or not len(right):
        raise SystemExit("The SV has no usable SNVs on one side; widen --flank")

    if args.genetic_map:
        rmap, rsrc = read_genetic_map(args.genetic_map, args.chrom, region.seq_len,
                                      args.recombination_rate)
    else:
        rmap, rsrc = args.recombination_rate, f"constant:{args.recombination_rate}"
    stats["recombination_source"] = rsrc

    log.info("tsinfer on %d inference + %d singleton sites", stats["inference_sites"],
             stats["singleton_sites"])
    t1 = time.time()
    ts = tsinfer.infer(sd, recombination_rate=rmap, mismatch_ratio=args.mismatch_ratio,
                       num_threads=args.threads)
    stats["tsinfer_seconds"] = round(time.time() - t1, 1)
    if not np.array_equal(ts.samples(), np.arange(n_haps)):
        raise SystemExit("Unexpected sample node layout from tsinfer")

    gap = region.excl_hi - region.excl_lo
    ts = tsdate.preprocess_ts(ts, minimum_gap=max(1_000_000, 2 * gap + 1))
    t1 = time.time()
    ts = tsdate.date(ts, mutation_rate=args.mutation_rate,
                     singletons_phased=args.singletons_phased)
    stats["tsdate_seconds"] = round(time.time() - t1, 1)
    muts_per_site = np.bincount(ts.mutations_site, minlength=ts.num_sites)
    stats.update(num_trees=ts.num_trees, num_nodes=ts.num_nodes, num_edges=ts.num_edges,
                 num_sites=ts.num_sites, num_mutations=ts.num_mutations,
                 frac_sites_multiple_mutations=round(float((muts_per_site > 1).mean()), 4))

    # ---- SV on the tree at its own position
    sv_pos = region.sv_start
    if not (ts.first().interval.left <= sv_pos < ts.sequence_length):
        raise SystemExit("SV position lies outside the inferred tree sequence")
    tree = ts.at(sv_pos)
    called = geno >= 0

    geno, placed = place_unphased(tree, geno, unphased_inds, hap_start, args.sv_ancestral)
    anc_state, trans = fitch(tree, geno, args.sv_ancestral)
    other = "alt" if anc_state == 0 else "ref"
    cost_flipped = cost(tree, geno, other)
    der = 1 - anc_state
    owner, under, rows, main = describe(tree, geno, anc_state, trans)
    rows = classify(rows, trans, main, anc_state, args.max_error_clade_haps,
                    args.recurrent_min_haps)
    in_main = np.zeros(n_haps, dtype=bool)
    in_main[under[main]] = True
    main_node = trans[main].node
    lo, hi, mid = branch_age(ts, tree, main_node)
    _, node_mn, node_vr = node_age(ts, main_node)

    # ---- suspicious individuals
    sus = []
    base_cost = len(trans)
    bad_haps = {}
    for r in rows:
        if r["role"] == "main" or r["role"] in ("recurrent_candidate", "large_reversion"):
            continue
        k = r["k"]
        mine = owner == k
        want = der if r["state"] == der else anc_state
        for h in np.flatnonzero(mine & (geno == want)):
            bad_haps[int(h)] = (r["role"], k)
    test_inds = sorted({int(hap_owner[h]) for h in bad_haps})
    swap_gain = {}
    for i in test_inds[:args.max_phase_tests]:
        if ploidy[i] != 2:
            continue
        a, b = hap_start[i], hap_start[i] + 1
        if geno[a] < 0 or geno[b] < 0 or geno[a] == geno[b]:
            continue
        alt = geno.copy()
        alt[a], alt[b] = geno[b], geno[a]
        swap_gain[i] = base_cost - cost(tree, alt, args.sv_ancestral)
    fixed = geno.copy()
    for i, gain in swap_gain.items():
        if gain > 0:
            a, b = hap_start[i], hap_start[i] + 1
            trial = fixed.copy()
            trial[a], trial[b] = fixed[b], fixed[a]
            if cost(tree, trial, args.sv_ancestral) < cost(tree, fixed, args.sv_ancestral):
                fixed = trial
    cost_after_phase = cost(tree, fixed, args.sv_ancestral)

    def gt_str(i):
        return "|".join(str(int(geno[hap_start[i] + j])) if geno[hap_start[i] + j] >= 0 else "."
                        for j in range(ploidy[i]))

    def main_str(i):
        return "|".join(str(der if in_main[hap_start[i] + j] else anc_state)
                        for j in range(ploidy[i]))

    for h, (role, k) in sorted(bad_haps.items()):
        i = int(hap_owner[h])
        gain = swap_gain.get(i, 0)
        issue = ("phasing_error" if gain > 0 else
                 "possible_false_positive" if geno[h] == 1 and der == 1 else
                 "possible_false_negative" if geno[h] == 0 and der == 1 else
                 "possible_false_negative" if geno[h] == 1 else "possible_false_positive")
        if role.startswith("ambiguous"):
            issue += "_weak"
        sus.append({"sample": samples[i], "haplotype": int(h - hap_start[i]),
                    "issue": issue, "clade_role": role, "clade_size_haps": rows[k]["haps_exclusive"],
                    "called_gt": gt_str(i), "main_clade_gt": main_str(i),
                    "other_hap_in_main_clade": bool(ploidy[i] == 2 and
                                                    in_main[hap_start[i] + (1 - (h - hap_start[i]))]),
                    "cost_gain_if_phase_swapped": int(gain)})
    for i, hap in placed:
        sus.append({"sample": samples[i], "haplotype": int(hap - hap_start[i]),
                    "issue": "unphased_call_placed_by_tree", "clade_role": "",
                    "clade_size_haps": np.nan, "called_gt": gt_str(i), "main_clade_gt": main_str(i),
                    "other_hap_in_main_clade": np.nan, "cost_gain_if_phase_swapped": np.nan})
    for i in np.flatnonzero([(geno[hap_start[i]:hap_start[i] + ploidy[i]] < 0).all()
                             for i in range(len(samples))]):
        if in_main[hap_start[i]:hap_start[i] + ploidy[i]].any():
            sus.append({"sample": samples[i], "haplotype": np.nan, "issue": "missing_call_in_main_clade",
                        "clade_role": "", "clade_size_haps": np.nan, "called_gt": gt_str(i),
                        "main_clade_gt": main_str(i), "other_hap_in_main_clade": np.nan,
                        "cost_gain_if_phase_swapped": np.nan})

    # ---- tree scan around the SV
    scan = []
    for t in ts.trees():
        if t.interval.right <= region.sv_start - args.scan_bp:
            continue
        if t.interval.left > region.sv_end + args.scan_bp:
            break
        a_, tr = fitch(t, geno, args.sv_ancestral)
        _, und, rr, mk = describe(t, geno, a_, tr)
        l_, h_, m_ = branch_age(ts, t, tr[mk].node)
        d_ = 1 - a_
        scan.append({"tree_index": t.index, "left": t.interval.left, "right": t.interval.right,
                     "contains_sv": t.interval.left <= sv_pos < t.interval.right,
                     "parsimony_cost": len(tr), "main_node": tr[mk].node,
                     "main_clade_derived_haps": rr[mk]["derived_called"],
                     "carriers_in_main_frac": rr[mk]["derived_called"] / max(1, int((geno == d_).sum())),
                     "age_gen_lower": l_, "age_gen_upper": h_, "age_gen_mid": m_})
    scan = pd.DataFrame(scan)
    # trees far from the SV carry a different genealogy; summarise only those that fit it well
    best = scan[scan.parsimony_cost <= scan.parsimony_cost.min() + 1]

    # ---- background SNV cost for comparable frequencies
    if len(holdout):
        bg, bg_split = background_costs(ts, held, args.recurrent_min_haps)
    else:
        bg, bg_split = np.array([]), np.array([], dtype=bool)
    pct = (float(((bg < base_cost).sum() + 0.5 * (bg == base_cost).sum()) / len(bg))
           if len(bg) else np.nan)
    split_rate = float(bg_split.mean()) if len(bg_split) else np.nan

    # ---- flags
    n_fp_fn = sum(1 for s in sus if s["issue"].startswith("possible_"))
    n_phase = len({s["sample"] for s in sus if s["issue"] == "phasing_error"})
    rec = [r for r in rows if r.get("role") == "recurrent_candidate"]
    rev = [r for r in rows if r.get("role") == "large_reversion"]
    low_density = (stats["sites_within_10kb_left"] < args.min_flank_sites or
                   stats["sites_within_10kb_right"] < args.min_flank_sites)
    within_bg = bool(len(bg) and pct < 0.9)
    issues = []
    if n_fp_fn:
        issues.append("GENOTYPE_ERROR")
    if n_phase:
        issues.append("PHASING_ERROR")
    if rec:
        common = np.isfinite(split_rate) and split_rate >= 0.05
        issues.append("POSSIBLE_RECURRENT" + ("(clade_splits_common_in_background)" if common else ""))
    if rev:
        issues.append("LARGE_REVERSION")
    on_root = tree.parent(main_node) == tskit.NULL
    polarity_suspect = args.sv_ancestral != "auto" and (cost_flipped < base_cost or on_root)
    if polarity_suspect:
        issues.append("POLARITY_SUSPECT")
    if low_density:
        issues.append("LOW_SITE_DENSITY")
    if anc.source == "reference":
        issues.append("REF_AS_ANCESTRAL")
    status = "PASS" if base_cost == 1 and not low_density else ";".join(issues) or "PASS"
    if base_cost > 1 and within_bg:
        status += "(within_background)"

    years = args.generation_time
    carriers_der = int((geno == der).sum())
    summary = {
        "sv_label": args.prefix, "sv_ids": ",".join(sv_info["sv_ids"]),
        "found_in_phased": ",".join(sv_info["found_phased"]),
        "found_in_full": ",".join(sv_info["found_full"]),
        "not_found": ",".join(sv_info["not_found"]),
        "chrom": args.chrom, "sv_start": region.sv_start, "sv_end": region.sv_end,
        "map_position": sv_pos, "sv_type": args.sv_type,
        "sv_ancestral_mode": args.sv_ancestral,
        "ancestral_state": "REF" if anc_state == 0 else "ALT",
        "dated_allele": "ALT" if der == 1 else "REF",
        "n_haplotypes": n_haps, "n_called_haplotypes": int(called.sum()),
        "derived_haplotypes": carriers_der,
        "derived_freq": carriers_der / max(1, int(called.sum())),
        "n_unphased_calls_placed": len(placed),
        "n_individuals_missing_sv_call": int(sum(1 for i in range(len(samples))
                                                 if (geno[hap_start[i]:hap_start[i] + ploidy[i]] < 0).all())),
        "tree_index": tree.index, "tree_left": tree.interval.left, "tree_right": tree.interval.right,
        "parsimony_cost": base_cost,
        "n_forward_mutations": sum(1 for r in rows if r["state"] == der),
        "n_back_mutations": sum(1 for r in rows if r["state"] != der),
        "cost_if_polarity_flipped": cost_flipped,
        "main_node": main_node, "main_clade_haps": rows[main]["haps_total"],
        "main_clade_derived_haps": rows[main]["derived_called"],
        "carriers_in_main_frac": rows[main]["derived_called"] / max(1, carriers_der),
        "main_clade_purity": rows[main]["derived_called"] /
                             max(1, rows[main]["derived_called"] + rows[main]["ancestral_called"]),
        "age_gen_lower": lo, "age_gen_upper": hi, "age_gen_mid": mid,
        "age_years_lower": lo * years, "age_years_upper": hi * years, "age_years_mid": mid * years,
        "main_node_posterior_mean_gen": node_mn,
        "main_node_posterior_sd_gen": float(np.sqrt(node_vr)) if np.isfinite(node_vr) else np.nan,
        "scan_n_trees": len(scan),
        "scan_min_cost": int(scan.parsimony_cost.min()),
        "scan_frac_trees_cost1": float((scan.parsimony_cost == 1).mean()),
        "scan_n_best_trees": len(best),
        "scan_best_age_gen_median": float(best.age_gen_mid.median()),
        "scan_best_age_gen_min": float(best.age_gen_mid.min()),
        "scan_best_age_gen_max": float(best.age_gen_mid.max()),
        "scan_frac_trees_same_main_node": float((scan.main_node == main_node).mean()),
        "background_n_sites": len(bg), "background_derived_count_range": f"{bg_lo:g}-{bg_hi:g}",
        "background_mean_cost": float(bg.mean()) if len(bg) else np.nan,
        "background_frac_cost_gt1": float((bg > 1).mean()) if len(bg) else np.nan,
        "background_frac_split_clades": split_rate,
        "cost_percentile_vs_background": pct,
        "n_suspect_genotype_haps": n_fp_fn, "n_phase_error_individuals": n_phase,
        "cost_after_phase_fixes": cost_after_phase,
        "n_recurrent_clades": len(rec),
        "recurrent_clade_carrier_haps": ",".join(str(r["derived_called"]) for r in rec),
        "n_large_reversions": len(rev),
        "flag_genotype_error": bool(n_fp_fn), "flag_phasing_error": bool(n_phase),
        "flag_recurrent": bool(rec), "flag_low_site_density": low_density,
        "flag_reference_as_ancestral": anc.source == "reference",
        "flag_singleton_sv": carriers_der == 1,
        "flag_polarity_suspect": polarity_suspect,
        "main_mutation_on_root": on_root,
        "inconsistency_within_background": within_bg if base_cost > 1 else np.nan,
        "qc_status": status,
    }

    # ---- write outputs
    p = args.prefix
    pd.DataFrame([summary]).to_csv(f"{p}.sv_age.tsv", sep="\t", index=False, na_rep="NA")
    mut = pd.DataFrame(rows)
    mut["parent_node"] = mut.node.map(tree.parent)
    mut["parent_mutation_node"] = mut.parent_mut.map(
        lambda k: trans[k].node if k != tskit.NULL else tskit.NULL)
    mut["age_gen_lower"] = mut.node.map(lambda n: branch_age(ts, tree, n)[0])
    mut["age_gen_upper"] = mut.node.map(lambda n: branch_age(ts, tree, n)[1])
    mut["age_gen_mid"] = mut.node.map(lambda n: branch_age(ts, tree, n)[2])
    mut["derived_state"] = np.where(mut.state == 1, "ALT", "REF")
    mut.drop(columns=["state", "k", "parent_mut"]).to_csv(f"{p}.sv_mutations.tsv", sep="\t", index=False,
                                            na_rep="NA")
    cols = ["sample", "haplotype", "issue", "clade_role", "clade_size_haps", "called_gt",
            "main_clade_gt", "other_hap_in_main_clade", "cost_gain_if_phase_swapped"]
    pd.DataFrame(sus, columns=cols).to_csv(f"{p}.suspicious_samples.tsv", sep="\t", index=False,
                                           na_rep="NA")
    pd.DataFrame({"sample": samples,
                  "called_gt": [gt_str(i) for i in range(len(samples))],
                  "main_clade_gt": [main_str(i) for i in range(len(samples))]}
                 ).to_csv(f"{p}.sv_genotypes.tsv.gz", sep="\t", index=False)
    scan.to_csv(f"{p}.tree_scan.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(sv_info["entries"]).to_csv(f"{p}.sv_entries.tsv", sep="\t", index=False)

    # ---- the ARG, with the SV added as a site
    tables = ts.dump_tables()
    site = tables.sites.add_row(position=sv_pos, ancestral_state="REF" if anc_state == 0 else "ALT")
    for r in rows:
        lo_, hi_, mid_ = branch_age(ts, tree, r["node"])
        t_mut = mid_ if np.isfinite(mid_) else lo_
        tables.mutations.add_row(site=site, node=r["node"], time=t_mut,
                                 derived_state="ALT" if r["state"] == 1 else "REF")
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tables.provenances.add_row(record=json.dumps({
        "software": {"name": "sv_age_arg"},
        "parameters": {k: v for k, v in vars(args).items()},
        "sv": {"site_id": int(site), "position": sv_pos, "ids": sv_info["sv_ids"],
               "chrom": args.chrom, "start": region.sv_start, "end": region.sv_end}}, default=str))
    out_ts = tables.tree_sequence()
    tszip.compress(out_ts, f"{p}.arg.trees.tsz")

    stats["total_seconds"] = round(time.time() - t0, 1)
    pd.DataFrame(list(stats.items()), columns=["metric", "value"]).to_csv(
        f"{p}.region_qc.tsv", sep="\t", index=False)
    log.info("SV %s: cost %d, age %.0f-%.0f generations (mid %.0f), status %s",
             args.prefix, base_cost, lo, hi, mid, status)
    for f in (f"{p}.samples", f"{p}.infer.samples"):
        if os.path.exists(f):
            os.remove(f)


if __name__ == "__main__":
    main()
