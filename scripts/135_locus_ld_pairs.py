#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""135_locus_ld_pairs.py -- pairwise r2 among gene lead SNPs and MR instruments,
computed from the saved standardised dosage matrix (1000G EUR, 503 samples).

Writes results/Coloc/01_tables/
  T10_gene_lead_r2.csv          r2 between every pair of gene lead SNPs
  T11_instrument_vs_lead_r2.csv r2 between every instrument and every gene lead
  T12_instrument_r2_matrix.csv  r2 among the six instruments
  T13_lead_blocks.csv           SNPs with r2 = 1.00 to the AQP3 lead (LD resolution)
ASCII only.
"""
import csv
import os

import numpy as np

ROOT = r"E:/another data file/GEObulk analysis project"
LD = os.path.join(ROOT, "results", "Coloc", "05_ld")
DAT = os.path.join(ROOT, "results", "Coloc", "01_data")
TAB = os.path.join(ROOT, "results", "Coloc", "01_tables")

GENES = ["AQP3", "RP11-176F3.7", "GALT", "DCTN3", "NOL6", "RP11-384P7.7", "BAG1"]
INSTR = ["rs10971815", "rs116972937", "rs12551364", "rs3814511", "rs573538", "rs626069"]


def main():
    Z = np.load(os.path.join(LD, "dosage_std.npy"))
    snps = [r["snp"] for r in csv.DictReader(open(os.path.join(LD, "dosage_snps.csv"),
                                                  newline="", encoding="utf-8"))]
    idx = {s: i for i, s in enumerate(snps)}
    R = (Z @ Z.T) / Z.shape[1]          # 4813 x 4813, fine for this size

    def r2(a, b):
        if a not in idx or b not in idx:
            return None
        return float(R[idx[a], idx[b]] ** 2)

    # gene leads: max |z1| within that gene's own SNP set (FinnGen exposure table)
    leads = {}
    for g in GENES:
        f = os.path.join(DAT, "coloc_%s_FinnGen.csv" % g.replace(".", "_"))
        if not os.path.exists(f):
            continue
        best, bz = None, 0.0
        for row in csv.DictReader(open(f, newline="", encoding="utf-8")):
            z = abs(float(row["z1"]))
            if z > bz:
                bz, best = z, row["snp"]
        leads[g] = (best, bz, float(row["pos"]))

    rows = []
    for a in GENES:
        for b in GENES:
            if a >= b:
                continue
            v = r2(leads[a][0], leads[b][0])
            rows.append({"gene_a": a, "snp_a": leads[a][0], "gene_b": b, "snp_b": leads[b][0],
                         "r2": "" if v is None else round(v, 4)})
    with open(os.path.join(TAB, "T10_gene_lead_r2.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["gene_a", "snp_a", "gene_b", "snp_b", "r2"])
        w.writeheader(); w.writerows(rows)
    print("gene lead SNPs:")
    for g, (s, z, p) in leads.items():
        print("  %-14s %-12s z=%6.2f pos=%.0f" % (g, s, z, p))
    print("\nr2 between gene leads:")
    for r in sorted(rows, key=lambda x: -(x["r2"] if x["r2"] != "" else -1)):
        print("  %-14s x %-14s %s" % (r["gene_a"], r["gene_b"], r["r2"]))

    rows = []
    for g in GENES:
        for ins in INSTR:
            v = r2(ins, leads[g][0])
            rows.append({"instrument": ins, "gene": g, "lead": leads[g][0],
                         "r2": "" if v is None else round(v, 4)})
    with open(os.path.join(TAB, "T11_instrument_vs_lead_r2.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["instrument", "gene", "lead", "r2"])
        w.writeheader(); w.writerows(rows)
    print("\nr2 instrument vs gene lead:")
    for g in GENES:
        vals = [(r["instrument"], r["r2"]) for r in rows if r["gene"] == g and r["r2"] != ""]
        print("  %-14s %s" % (g, "  ".join("%s=%.3f" % (a, b) for a, b in vals)))

    M = np.zeros((len(INSTR), len(INSTR)))
    for i, a in enumerate(INSTR):
        for j, b in enumerate(INSTR):
            v = r2(a, b)
            M[i, j] = np.nan if v is None else v
    with open(os.path.join(TAB, "T12_instrument_r2_matrix.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([""] + INSTR)
        for i, a in enumerate(INSTR):
            w.writerow([a] + ["%.3f" % x if not np.isnan(x) else "NA" for x in M[i]])
    print("\nr2 among instruments (upper triangle):")
    for i, a in enumerate(INSTR):
        print("  %-12s %s" % (a, "  ".join(
            "%s=%.3f" % (INSTR[j], M[i, j]) for j in range(i + 1, len(INSTR)))))

    # LD block containing the AQP3 lead
    lp = leads["AQP3"][0]
    blk = [(s, R[idx[lp], idx[s]]) for s in snps if abs(R[idx[lp], idx[s]]) > 0.999]
    with open(os.path.join(TAB, "T13_lead_blocks.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene", "lead", "snp", "r"])
        for s, r in sorted(blk, key=lambda x: -x[1]):
            w.writerow(["AQP3", lp, s, round(r, 4)])
    print("\nAQP3 lead LD block (|r| > 0.999): %d SNPs" % len(blk))
    for s, r in sorted(blk, key=lambda x: -x[1]):
        print("  %-12s r=%.3f" % (s, r))


if __name__ == "__main__":
    main()
