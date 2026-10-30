#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""126_locus_coloc_harmonise.py -- build per-gene x per-outcome coloc input tables.

For every gene of the 9p13.3 block and every outcome (HERMES, FinnGen) writes
  results/Coloc/01_data/coloc_<gene>_<outcome>.csv
with columns: snp, pos, beta1, varbeta1, beta2, varbeta2, z1, z2, maf

Exposure betas are reconstructed from eQTLGen Z scores:
  se = 1 / sqrt(2*MAF*(1-MAF)*(N + Z^2));  beta = Z * se
(Giambartolomei 2014 / Zhu 2016 convention for z-only eQTL releases).

Note on palindromic variants: the approximate Bayes factor depends on beta^2,
so allele flips do not change it; palindromic SNPs are therefore RETAINED. Only
variants whose alleles cannot be reconciled with the eQTL alleles are dropped.
ASCII only.
"""
import csv
import math
import os
import sys
from collections import defaultdict

ROOT = r"E:/another data file/GEObulk analysis project"
DAT = os.path.join(ROOT, "results", "Coloc", "01_data")
GENES = ["AQP3", "RP11-176F3.7", "DCTN3", "GALT", "NUDT2", "RP11-384P7.7", "BAG1", "NOL6"]
OUTCOMES = {
    "HERMES": {"file": "hermes_region.csv", "n": 977323, "s": 47309 / 977323.0},
    "FinnGen": {"file": "finngen_region.csv", "n": 453733, "s": 33250 / 453733.0},
}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def read_csv(name):
    with open(os.path.join(DAT, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    eqtl = read_csv("eqtl_region_long.csv")
    maf = {r["rsid"]: r for r in read_csv("maf_region.csv")}

    by_gene = defaultdict(list)
    for r in eqtl:
        by_gene[r["gene"]].append(r)

    summary = []
    for oname, cfg in OUTCOMES.items():
        out_rows = {}
        for r in read_csv(cfg["file"]):
            out_rows.setdefault(r["snp"], r)  # first record wins (no duplicate rsIDs expected)
        for gene in GENES:
            rows = by_gene.get(gene, [])
            recs, n_no_maf, n_no_out, n_bad_allele, n_palindromic = [], 0, 0, 0, 0
            for r in rows:
                snp = r["snp"]
                m = maf.get(snp)
                if not m or m["maf"] in ("", "NA"):
                    n_no_maf += 1
                    continue
                o = out_rows.get(snp)
                if o is None:
                    n_no_out += 1
                    continue
                a1, a2 = r["a1"].upper(), r["a2"].upper()
                ea, oa = o["ea"].upper(), o["oa"].upper()
                pal = (a1 in COMP and COMP[a1] == a2)
                if ea == a1 and oa == a2:
                    flip = 1.0
                elif ea == a2 and oa == a1:
                    flip = -1.0
                else:
                    n_bad_allele += 1
                    continue
                if pal:
                    n_palindromic += 1
                try:
                    z = float(r["z"])
                    n1 = float(r["n"])
                    mf = float(m["maf"])
                    b2 = float(o["beta"])
                    se2 = float(o["se"])
                except ValueError:
                    n_bad_allele += 1
                    continue
                if not (0 < mf < 1) or se2 <= 0:
                    n_bad_allele += 1
                    continue
                se1 = 1.0 / math.sqrt(2 * mf * (1 - mf) * (n1 + z * z))
                recs.append({
                    "snp": snp, "pos": int(r["pos"]),
                    "beta1": z * se1, "varbeta1": se1 ** 2,
                    "beta2": flip * b2, "varbeta2": se2 ** 2,
                    "z1": z, "z2": b2 / se2, "maf": mf,
                })
            recs.sort(key=lambda x: x["pos"])
            path = os.path.join(DAT, "coloc_%s_%s.csv" % (gene.replace(".", "_"), oname))
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=["snp", "pos", "beta1", "varbeta1",
                                                  "beta2", "varbeta2", "z1", "z2", "maf"])
                w.writeheader()
                w.writerows(recs)
            # lead SNPs
            lead1 = max(recs, key=lambda x: abs(x["z1"])) if recs else None
            lead2 = max(recs, key=lambda x: abs(x["z2"])) if recs else None
            summary.append({
                "gene": gene, "outcome": oname, "n_snp": len(recs),
                "drop_no_maf": n_no_maf, "drop_no_outcome": n_no_out,
                "drop_allele": n_bad_allele, "n_palindromic": n_palindromic,
                "lead_eqtl_snp": lead1["snp"] if lead1 else "",
                "lead_eqtl_pos": lead1["pos"] if lead1 else "",
                "lead_eqtl_z": round(lead1["z1"], 2) if lead1 else "",
                "lead_out_snp": lead2["snp"] if lead2 else "",
                "lead_out_pos": lead2["pos"] if lead2 else "",
                "lead_out_z": round(lead2["z2"], 2) if lead2 else "",
                "lead_distance_bp": abs(lead1["pos"] - lead2["pos"]) if (lead1 and lead2) else "",
            })
            print("%-14s %-8s n=%5d (noMAF %4d, noOut %4d, badAllele %4d, palindromic %4d)"
                  % (gene, oname, len(recs), n_no_maf, n_no_out, n_bad_allele, n_palindromic), flush=True)

    out = os.path.join(DAT, "coloc_input_summary.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
