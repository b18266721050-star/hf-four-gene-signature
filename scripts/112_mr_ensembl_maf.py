#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""112_mr_ensembl_maf.py -- look up 1000G/gnomAD minor allele frequency for instrument SNPs
via the Ensembl REST variation endpoint (batch POST, 150 ids per call).
Output: results/MR/01_data/T07_snp_maf.csv  (rsid, maf, allele_string, n_alleles, minor_allele)
ASCII only.
"""
import csv, json, os, sys, time
import urllib.request as ul

ROOT = r"E:/another data file/GEObulk analysis project"
OUT = os.path.join(ROOT, "results", "MR", "01_data")
BASE = "https://rest.ensembl.org/variation/homo_sapiens"


def post(ids):
    body = json.dumps({"ids": ids}).encode()
    req = ul.Request(BASE, data=body, headers={
        "Content-Type": "application/json", "Accept": "application/json"})
    for k in range(5):
        try:
            with ul.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            sys.stderr.write("retry %d: %s\n" % (k, e))
            time.sleep(3 + 4 * k)
    return {}


def main():
    rs = []
    for fn in ["T02_eqtlgen_instruments_pruned.csv", "T05_gtex_instruments_pruned.csv"]:
        p = os.path.join(OUT, fn)
        if os.path.exists(p):
            with open(p, newline="", encoding="utf-8") as f:
                for row in csv.DictReader(f):
                    rs.append(row["rsid"])
    rs = sorted(set(rs))
    print("instrument rsIDs:", len(rs), flush=True)
    rows, miss = [], []
    for i in range(0, len(rs), 150):
        chunk = rs[i:i + 150]
        d = post(chunk)
        for r in chunk:
            v = d.get(r)
            if not v:
                miss.append(r)
                continue
            al = v.get("allele_string") or ((v.get("mappings") or [{}])[0].get("allele_string") or "")
            n_alt = len([a for a in al.split("/") if a.strip()]) if al else 0
            rows.append({"rsid": r, "maf": v.get("MAF", ""),
                         "allele_string": al, "n_alleles": n_alt,
                         "minor_allele": v.get("minor_allele", ""),
                         "ambiguity": v.get("ambiguity", "")})
        print("chunk %d ok, cumulative %d" % (i // 150 + 1, len(rows)), flush=True)
        time.sleep(0.6)
    out = os.path.join(OUT, "T07_snp_maf.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["rsid", "maf", "allele_string", "n_alleles",
                                          "minor_allele", "ambiguity"])
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, len(rows), "rows;", len(miss), "not found", flush=True)
    if miss:
        print("missing:", ",".join(miss[:20]), flush=True)


if __name__ == "__main__":
    main()
