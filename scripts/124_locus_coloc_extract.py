#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""124_locus_coloc_extract.py -- extract the 9p13.3 AQP3 locus for colocalisation.

Outputs (results/Coloc/01_data/):
  eqtl_region_long.csv   eQTLGen cis-eQTLs for the 8 genes of the block (all FDR-passing pairs)
  hermes_region.csv      HERMES HF associations on chr9 in the window
  finngen_region.csv     FinnGen R11 I9_HEARTFAIL associations on chr9 in the window

Region: chr9 32.5-35.5 Mb in BOTH builds (the 9p13.3 block does not move materially
between GRCh37 and GRCh38; matching is done on rsID, coordinates are only a prefilter).
ASCII only.
"""
import csv
import gzip
import io
import os
import sys
import time

ROOT = r"E:/another data file/GEObulk analysis project"
OUT = os.path.join(ROOT, "results", "Coloc", "01_data")
GENES = ["AQP3", "DCTN3", "GALT", "NUDT2", "RP11-176F3.7", "RP11-384P7.7", "BAG1", "NOL6"]
CHR = "9"
LO, HI = 32000000, 37000000


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def extract_eqtl():
    src = os.path.join(ROOT, "raw_data", "mr", "eqtlgen_cis_eqtl_FDR005.txt.gz")
    dst = os.path.join(OUT, "eqtl_region_long.csv")
    want = set(GENES)
    n_in = n_out = 0
    with gzip.open(src, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        idx = {k: i for i, k in enumerate(head)}
        i_snp, i_chr, i_pos = idx["SNP"], idx["SNPChr"], idx["SNPPos"]
        i_a1, i_a2, i_z = idx["AssessedAllele"], idx["OtherAllele"], idx["Zscore"]
        i_gs, i_n = idx["GeneSymbol"], idx["NrSamples"]
        with open(dst, "w", newline="", encoding="utf-8") as out:
            w = csv.writer(out)
            w.writerow(["gene", "snp", "chr", "pos", "a1", "a2", "z", "n"])
            for row in rd:
                n_in += 1
                if len(row) < len(head):
                    continue
                if row[i_gs] not in want:
                    continue
                w.writerow([row[i_gs], row[i_snp], row[i_chr], row[i_pos],
                            row[i_a1], row[i_a2], row[i_z], row[i_n]])
                n_out += 1
    log("eqtlgen: %d rows scanned, %d kept -> %s" % (n_in, n_out, dst))


def extract_hermes():
    src = os.path.join(ROOT, "raw_data", "mr", "hermes_HF_GCST009541_h.tsv.gz")
    dst = os.path.join(OUT, "hermes_region.csv")
    n_in = n_out = 0
    with gzip.open(src, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        idx = {k: i for i, k in enumerate(head)}
        i_rs, i_chr = idx["hm_rsid"], idx["hm_chrom"]
        i_pos = idx["hm_pos"]
        i_ea, i_oa = idx["hm_effect_allele"], idx["hm_other_allele"]
        i_b, i_se = idx["hm_beta"], idx["standard_error"]
        i_braw = idx["beta"]
        i_p = idx["p_value"]
        # raw (unharmonised) columns, used when the hm_* fields are NA
        i_rchr, i_rpos = idx["chromosome"], idx["base_pair_location"]
        i_rea, i_roa = idx["effect_allele"], idx["other_allele"]
        i_rb = idx["beta"]

        def pick(row, i_h, i_r):
            v = row[i_h]
            return v if v not in ("", "NA", "nan") else row[i_r]

        with open(dst, "w", newline="", encoding="utf-8") as out:
            w = csv.writer(out)
            w.writerow(["snp", "chr", "pos", "ea", "oa", "beta", "se", "p"])
            for row in rd:
                n_in += 1
                if len(row) < len(head):
                    continue
                chrv = pick(row, i_chr, i_rchr)
                if chrv != CHR:
                    continue
                try:
                    pos = int(float(pick(row, i_pos, i_rpos)))
                except ValueError:
                    continue
                if not (LO <= pos <= HI):
                    continue
                if not row[i_rs].startswith("rs"):
                    continue
                ea = pick(row, i_ea, i_rea)
                oa = pick(row, i_oa, i_roa)
                b = pick(row, i_b, i_rb)
                # guard: harmonised beta must equal the raw beta, otherwise the
                # standard_error column refers to a different allele coding
                try:
                    if abs(float(b) - float(row[i_braw])) > 1e-6:
                        continue
                    pos = int(pos)
                except ValueError:
                    continue
                w.writerow([row[i_rs], chrv, pos, ea, oa, b, row[i_se], row[i_p]])
                n_out += 1
    log("hermes: %d rows scanned, %d kept -> %s" % (n_in, n_out, dst))


def extract_finngen():
    src = os.path.join(ROOT, "raw_data", "mr", "finngen_R11_I9_HEARTFAIL.gz")
    dst = os.path.join(OUT, "finngen_region.csv")
    n_in = n_out = 0
    with gzip.open(src, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(head)}
        i_chr, i_pos = idx["chrom"], idx["pos"]
        i_ref, i_alt, i_rs = idx["ref"], idx["alt"], idx["rsids"]
        i_b, i_se, i_p = idx["beta"], idx["sebeta"], idx["pval"]
        i_af = idx["af_alt"]
        out_f = open(dst, "w", newline="", encoding="utf-8")
        w = csv.writer(out_f)
        w.writerow(["snp", "chr", "pos", "ea", "oa", "beta", "se", "p", "af_alt"])
        for line in fh:
            n_in += 1
            row = line.rstrip("\n").split("\t")
            if len(row) < len(head):
                continue
            if row[i_chr] != CHR:
                continue
            try:
                pos = int(row[i_pos])
            except ValueError:
                continue
            if not (LO <= pos <= HI):
                continue
            for rs in row[i_rs].split(","):
                rs = rs.strip()
                if rs.startswith("rs"):
                    w.writerow([rs, row[i_chr], row[i_pos], row[i_alt], row[i_ref],
                                row[i_b], row[i_se], row[i_p], row[i_af]])
                    n_out += 1
                    break
            if n_in % 2000000 == 0:
                log("finngen scanned %d" % n_in)
        out_f.close()
    log("finngen: %d rows scanned, %d kept -> %s" % (n_in, n_out, dst))


def main():
    os.makedirs(OUT, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "eqtl"):
        extract_eqtl()
    if which in ("all", "hermes"):
        extract_hermes()
    if which in ("all", "finngen"):
        extract_finngen()


if __name__ == "__main__":
    main()
