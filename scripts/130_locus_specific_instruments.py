#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""130_locus_specific_instruments.py -- repeat the AQP3 cis-MR with gene-specific
instruments only.

Gene specificity is taken from the project's own locus check (T22_locus_cis_targets.csv,
p < 5e-8 for a cis target): an instrument is "AQP3-specific" when no other gene of the
9p13.3 block reaches that threshold for the same SNP.

HERMES betas/SEs are read from the pre-harmonised instrument table (T10) so that the
all-instrument estimate stays identical to the published MR; FinnGen betas/SEs are
taken from the FinnGen R11 region file with allele harmonisation (this reproduces the
published OR 0.787 exactly, which is the check that the two paths agree).

Writes results/Coloc/01_tables/T04_instrument_specificity.csv
       results/Coloc/01_tables/T05_instrument_gene_matrix.csv
ASCII only.
"""
import math
import os

import pandas as pd

ROOT = r"E:/another data file/GEObulk analysis project"
DAT = os.path.join(ROOT, "results", "Coloc", "01_data")
MRD = os.path.join(ROOT, "results", "MR")
TAB = os.path.join(ROOT, "results", "Coloc", "01_tables")
OUTS = {"HERMES": ("hermes_region.csv", 47309, 977323),
        "FinnGen": ("finngen_region.csv", 33250, 453733)}
BLOCK = {"AQP3", "DCTN3", "GALT", "NUDT2", "RP11-176F3.7", "RP11-384P7.7", "BAG1", "NOL6"}


def ivw(bx, by, sy):
    w = [(bx[i] ** 2) / (sy[i] ** 2) for i in range(len(bx))]
    b = sum(w[i] * (by[i] / bx[i]) for i in range(len(bx))) / sum(w)
    se = math.sqrt(1.0 / sum(w))
    z = b / se
    p = math.erfc(abs(z) / math.sqrt(2))
    q = sum(w[i] * ((by[i] / bx[i]) - b) ** 2 for i in range(len(bx)))
    return b, se, p, q, len(bx)


def wald(bx, by, sy):
    b = by / bx
    se = sy / abs(bx)
    p = math.erfc(abs(b / se) / math.sqrt(2))
    return b, se, p, 0.0, 1


def main():
    inst = pd.read_csv(os.path.join(MRD, "01_data", "T10_harmonised_snps.csv"))
    inst = inst[(inst.exposure == "eQTLGen_blood") & (inst.gene == "AQP3")].reset_index(drop=True)
    t22 = pd.read_csv(os.path.join(MRD, "01_tables", "T22_locus_cis_targets.csv"))
    t22 = t22[(t22.gene_assayed == "AQP3") & (t22.p_target < 5e-8)]

    mat = []
    for rs in inst.rsid:
        others = sorted(set(t22[(t22.rsid == rs) & (t22.target_gene != "AQP3")].target_gene))
        mat.append({"rsid": rs, "n_other_genes_p5e8": len(others),
                    "other_genes": ";".join(others),
                    "specific": len(others) == 0})
    mat = pd.DataFrame(mat)
    mat = mat.merge(inst[["rsid", "fstat"]], on="rsid", how="left")
    mat.to_csv(os.path.join(TAB, "T05_instrument_gene_matrix.csv"), index=False)
    print(mat.to_string(index=False))
    specific = set(mat.loc[mat.specific, "rsid"])
    print("\nAQP3-specific:", sorted(specific))

    rows = []
    for oname, (fn, ncase, ntot) in OUTS.items():
        bx, by, sy, used = [], [], [], []
        if oname == "HERMES":
            for _, r in inst.iterrows():
                bx.append(float(r.bx)); by.append(float(r.by)); sy.append(float(r.sy)); used.append(r.rsid)
        else:
            reg = pd.read_csv(os.path.join(DAT, fn)).drop_duplicates("snp").set_index("snp")
            for _, r in inst.iterrows():
                if r.rsid not in reg.index:
                    continue
                o = reg.loc[r.rsid]
                ea, oa = str(o.ea).upper(), str(o.oa).upper()
                a1, a2 = str(r.effect_allele).upper(), str(r.other_allele).upper()
                if ea == a1 and oa == a2:
                    flip = 1.0
                elif ea == a2 and oa == a1:
                    flip = -1.0
                else:
                    continue
                bx.append(float(r.bx)); by.append(flip * float(o.beta)); sy.append(float(o.se)); used.append(r.rsid)

        def add(label, idx):
            bxs = [bx[i] for i in idx]; bys = [by[i] for i in idx]; ss = [sy[i] for i in idx]
            if not idx:
                rows.append({"outcome": oname, "set": label, "n_snp": 0, "snps": "",
                             "OR": "", "lo": "", "hi": "", "p": "", "Q": ""})
                return
            if len(idx) == 1:
                b, se, p, q, n = wald(bxs[0], bys[0], ss[0]); tag = label + " (wald)"
            else:
                b, se, p, q, n = ivw(bxs, bys, ss); tag = label
            rows.append({"outcome": oname, "set": tag, "n_snp": n,
                         "snps": ";".join(used[i] for i in idx),
                         "OR": round(math.exp(b), 4), "lo": round(math.exp(b - 1.96 * se), 4),
                         "hi": round(math.exp(b + 1.96 * se), 4), "p": p, "Q": round(q, 3)})

        add("all", list(range(len(used))))
        add("AQP3-specific only", [i for i, r in enumerate(used) if r in specific])
        add("shared only", [i for i, r in enumerate(used) if r not in specific])

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TAB, "T04_instrument_specificity.csv"), index=False)
    print("\n" + out.to_string(index=False))


if __name__ == "__main__":
    main()
