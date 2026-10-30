#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""131_ld_from_vcf.py -- genotypes for the 9p13.3 SNP set from the 1000G phase 3
chr9 VCF (EUR only) and LD (r) matrices per gene window.

Matching is by POSITION + ALLELES, not by rsID: the phase 3 VCF carries "." in the
ID column for most records. eQTLGen and 1000G phase 3 are both GRCh37, so positions
are directly comparable.

Outputs (results/Coloc/05_ld/)
  ld_<gene>.csv, ld_<gene>_snps.csv   LD r matrix and SNP order (ALT orientation)
  align.csv                           snp, ref, alt, eqtl_a1, flip (+1/-1)
  ld_summary.csv                      SNP counts before/after matching
ASCII only.
"""
import csv
import gzip
import os
import sys

import numpy as np

ROOT = r"E:/another data file/GEObulk analysis project"
LD = os.path.join(ROOT, "raw_data", "ld")
DAT = os.path.join(ROOT, "results", "Coloc", "01_data")
OUT = os.path.join(ROOT, "results", "Coloc", "05_ld")
VCF = os.path.join(LD, "all_chr9.vcf.gz")
PANEL = os.path.join(LD, "panel.txt")
LO, HI = 32000000, 37000000
MAX_MISS = 0.10
MAX_SNP_LD = 2000

NL = b"\n"
TAB = b"\t"


def eur_samples():
    keep = set()
    with open(PANEL, encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["super_pop"] == "EUR":
                keep.add(r["sample"])
    return keep


def wanted_by_pos():
    """pos -> list of [snp, a1, a2, set(genes)]"""
    d = {}
    with open(os.path.join(DAT, "eqtl_region_long.csv"), newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pos = int(row["pos"])
            a1, a2 = row["a1"].upper(), row["a2"].upper()
            if len(a1) != 1 or len(a2) != 1:
                continue
            hit = None
            for rec in d.get(pos, []):
                if rec[0] == row["snp"]:
                    hit = rec
                    break
            if hit is None:
                rec = [row["snp"], a1, a2, set([row["gene"]])]
                d.setdefault(pos, []).append(rec)
            else:
                hit[3].add(row["gene"])
    return d


def main():
    os.makedirs(OUT, exist_ok=True)
    eur = eur_samples()
    want = wanted_by_pos()
    n_target = sum(len(v) for v in want.values())
    print("EUR samples: %d, target SNPs: %d at %d positions" % (len(eur), n_target, len(want)),
          flush=True)

    snps, mat, alleles, genes_of = [], [], {}, {}
    state = {"win": 0, "hit": 0, "stop": 0}
    cols = None

    def handle(line):
        if not line.startswith(b"9" + TAB):
            return
        f = line.split(TAB, 6)
        if len(f) < 7:
            return
        try:
            pos = int(f[1])
        except ValueError:
            return
        if pos < LO:
            return
        if pos > HI:
            state["stop"] = 1
            return
        state["win"] += 1
        cand = want.get(pos)
        if not cand:
            return
        ref, alt = f[3].decode().upper(), f[4].decode().upper()
        if len(ref) != 1 or len(alt) != 1:
            return
        rec = None
        for r in cand:
            if (r[1] == alt and r[2] == ref) or (r[1] == ref and r[2] == alt):
                rec = r
                break
        if rec is None:
            return
        gt = f[6].split(TAB)
        if len(gt) < len(cols):
            return
        d = np.empty(len(cols), dtype=np.float32)
        miss = 0
        for k, ci in enumerate(cols):
            g = gt[ci]
            if g == b"0|0":
                d[k] = 0.0
            elif g == b"0|1" or g == b"1|0":
                d[k] = 1.0
            elif g == b"1|1":
                d[k] = 2.0
            else:
                d[k] = np.nan
                miss += 1
        if miss > MAX_MISS * len(cols):
            return
        snps.append(rec[0])
        mat.append(d)
        alleles[rec[0]] = (ref, alt)
        genes_of[rec[0]] = rec[3]
        state["hit"] += 1
        if state["hit"] % 500 == 0:
            print("  matched %d" % state["hit"], flush=True)

    with gzip.open(VCF, "rb") as fh:
        for line in fh:
            if line.startswith(b"##"):
                continue
            if line.startswith(b"#CHROM"):
                names = line.decode().rstrip().split("\t")[9:]
                cols = [i for i, s in enumerate(names) if s in eur]
                print("header samples: %d, EUR columns: %d" % (len(names), len(cols)), flush=True)
                break
        if cols is None:
            sys.exit("no #CHROM header found")

        chunk = b""
        while not state["stop"]:
            buf = fh.read(1 << 25)
            if not buf:
                break
            buf = chunk + buf
            lines = buf.split(NL)
            chunk = lines.pop()
            for line in lines:
                handle(line)
                if state["stop"]:
                    break

    print("in window %d, matched %d" % (state["win"], state["hit"]), flush=True)
    if not snps:
        sys.exit("nothing matched - check window / alleles")

    M = np.vstack(mat)
    mu = np.nanmean(M, axis=1, keepdims=True)
    M = np.where(np.isnan(M), mu, M)

    rows = []
    for i, s in enumerate(snps):
        ref, alt = alleles[s]
        rows.append({"snp": s, "ref": ref, "alt": alt})
    # flip needs the eQTLGen effect allele; recover it from the original table
    a1map = {}
    with open(os.path.join(DAT, "eqtl_region_long.csv"), newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            a1map.setdefault(row["snp"], row["a1"].upper())
    keep = []
    for i, s in enumerate(snps):
        ref, alt = alleles[s]
        a1 = a1map.get(s)
        if a1 == alt:
            flip = 1
        elif a1 == ref:
            flip = -1
        else:
            continue
        rows[i]["eqtl_a1"] = a1
        rows[i]["flip"] = flip
        keep.append(i)
    print("allele-reconciled SNPs: %d / %d" % (len(keep), len(snps)), flush=True)

    M = M[keep]
    snps = [snps[i] for i in keep]
    rows = [rows[i] for i in keep]
    with open(os.path.join(OUT, "align.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["snp", "ref", "alt", "eqtl_a1", "flip"])
        w.writeheader()
        w.writerows(rows)

    Z = M - M.mean(axis=1, keepdims=True)
    sd = Z.std(axis=1)
    Z = Z / np.where(sd == 0, 1, sd)[:, None]

    np.save(os.path.join(OUT, "dosage_std.npy"), Z.astype(np.float32))
    with open(os.path.join(OUT, "dosage_snps.csv"), "w", newline="", encoding="utf-8") as fh:
        ww = csv.writer(fh)
        ww.writerow(["snp"])
        for s in snps:
            ww.writerow([s])

    idx = {s: i for i, s in enumerate(snps)}
    by_gene = {}
    for s in snps:
        for g in genes_of[s]:
            by_gene.setdefault(g, []).append(s)

    summary = []
    with open(os.path.join(OUT, "ld_summary.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "n_1000g", "ld_dim"])
        w.writeheader()
        for g, ss in sorted(by_gene.items()):
            if len(ss) < 20 or len(ss) > MAX_SNP_LD:
                w.writerow({"gene": g, "n_1000g": len(ss), "ld_dim": 0})
                print("  skip LD %-14s (%d SNPs)" % (g, len(ss)), flush=True)
                continue
            ii = [idx[s] for s in ss]
            Zg = Z[ii]
            R = (Zg @ Zg.T) / Zg.shape[1]
            np.fill_diagonal(R, 1.0)
            name = g.replace(".", "_")
            np.savetxt(os.path.join(OUT, "ld_%s.csv" % name), R, delimiter=",", fmt="%.5f")
            with open(os.path.join(OUT, "ld_%s_snps.csv" % name), "w", newline="", encoding="utf-8") as fh:
                ww = csv.writer(fh)
                ww.writerow(["snp"])
                for s in ss:
                    ww.writerow([s])
            w.writerow({"gene": g, "n_1000g": len(ss), "ld_dim": R.shape[0]})
            print("  LD %-14s %4d" % (g, R.shape[0]), flush=True)


if __name__ == "__main__":
    main()
