#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""125_locus_coloc_maf.py -- 1000G/gnomAD MAF for the 9p13.3 locus SNPs (Ensembl REST batch).

Reads  results/Coloc/01_data/eqtl_region_long.csv
Writes results/Coloc/01_data/maf_region.csv  (rsid, maf, allele_string, n_alleles)

Ensembl batch POST is limited to 200 ids per call; the allele_string lives in
mappings[0] for some records (see 112_mr_ensembl_maf.py). ASCII only.
"""
import csv
import json
import os
import sys
import time
import urllib.request as ul

ROOT = r"E:/another data file/GEObulk analysis project"
DAT = os.path.join(ROOT, "results", "Coloc", "01_data")
BASE = "https://rest.ensembl.org/variation/homo_sapiens"


def post(ids):
    body = json.dumps({"ids": ids}).encode()
    req = ul.Request(BASE, data=body, headers={
        "Content-Type": "application/json", "Accept": "application/json"})
    for k in range(5):
        try:
            with ul.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            sys.stderr.write("retry %d (%s): %s\n" % (k, len(ids), e))
            time.sleep(3 + 4 * k)
    return {}


def main():
    snps = set()
    with open(os.path.join(DAT, "eqtl_region_long.csv"), newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            snps.add(row["snp"])
    snps = sorted(snps)
    print("locus rsIDs:", len(snps), flush=True)

    rows, miss = [], []
    for i in range(0, len(snps), 150):
        chunk = snps[i:i + 150]
        d = post(chunk)
        for r in chunk:
            v = d.get(r)
            if not v:
                miss.append(r)
                continue
            al = v.get("allele_string") or ((v.get("mappings") or [{}])[0].get("allele_string") or "")
            n_alt = len([a for a in al.split("/") if a.strip()]) if al else 0
            rows.append({"rsid": r, "maf": v.get("MAF", ""), "allele_string": al,
                         "n_alleles": n_alt, "minor_allele": v.get("minor_allele", "")})
        if (i // 150 + 1) % 5 == 0:
            print("  batch %d, cumulative %d" % (i // 150 + 1, len(rows)), flush=True)
        time.sleep(0.5)

    out = os.path.join(DAT, "maf_region.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["rsid", "maf", "allele_string", "n_alleles", "minor_allele"])
        w.writeheader()
        w.writerows(rows)
    print("wrote %s: %d rows, %d not found" % (out, len(rows), len(miss)), flush=True)
    if miss:
        print("first missing:", ",".join(miss[:10]))


if __name__ == "__main__":
    main()
