#!/usr/bin/env python3
"""
annotate_snvs.py
----------------
Reads bcftools query output from stdin, adds computed fields,
and writes an enriched TSV to stdout.

Added columns:
  SAMPLE, COHORT, VARIANT_TYPE, SNV_CLASS (ts/tv/indel),
  AF_APPROX (from AD), HET_FLAG, PASS_FLAG, CHROM_TYPE
"""

import argparse
import sys

import pandas as pd


TRANSITION_PAIRS = {
    ("A", "G"), ("G", "A"),
    ("C", "T"), ("T", "C"),
}


def snv_class(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        pair = (ref.upper(), alt.upper())
        return "Ti" if pair in TRANSITION_PAIRS else "Tv"
    return "INDEL"


def approx_af(ad_field: str) -> float:
    """Compute allele frequency from AD field (ref,alt)."""
    try:
        parts = str(ad_field).split(",")
        if len(parts) < 2:
            return float("nan")
        ref_d, alt_d = int(parts[0]), int(parts[1])
        total = ref_d + alt_d
        return round(alt_d / total, 4) if total > 0 else float("nan")
    except (ValueError, ZeroDivisionError):
        return float("nan")


def chrom_type(chrom: str) -> str:
    c = str(chrom).lower()
    if c == "chrx":
        return "chrX"
    if c == "chry":
        return "chrY"
    if c == "chrm" or c == "chrmt":
        return "chrMT"
    return "autosome"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample",  required=True)
    p.add_argument("--cohort",  required=True)
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(
        sys.stdin,
        sep="\t",
        dtype=str,
    )

    expected = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "GT", "DP", "GQ", "AD"]
    for col in expected:
        if col not in df.columns:
            df[col] = "."

    df["SAMPLE"]       = args.sample
    df["COHORT"]       = args.cohort
    df["QUAL"]         = pd.to_numeric(df["QUAL"], errors="coerce")
    df["DP"]           = pd.to_numeric(df["DP"],   errors="coerce")
    df["GQ"]           = pd.to_numeric(df["GQ"],   errors="coerce")

    df["SNV_CLASS"]    = df.apply(lambda r: snv_class(str(r["REF"]), str(r["ALT"])), axis=1)
    df["AF_APPROX"]    = df["AD"].apply(approx_af)
    df["HET_FLAG"]     = df["GT"].apply(lambda g: 1 if "0/1" in str(g) or "0|1" in str(g) else 0)
    df["HOM_ALT_FLAG"] = df["GT"].apply(lambda g: 1 if str(g) in ("1/1", "1|1") else 0)
    df["PASS_FLAG"]    = df["FILTER"].apply(lambda f: 1 if str(f).upper() in ("PASS", ".") else 0)
    df["CHROM_TYPE"]   = df["CHROM"].apply(chrom_type)
    df["VARIANT_TYPE"] = "SNV"

    # Reorder
    out_cols = [
        "SAMPLE", "COHORT", "CHROM", "CHROM_TYPE", "POS",
        "ID", "REF", "ALT", "QUAL", "FILTER", "PASS_FLAG",
        "GT", "DP", "GQ", "AD", "AF_APPROX",
        "HET_FLAG", "HOM_ALT_FLAG", "SNV_CLASS", "VARIANT_TYPE",
    ]
    df = df[[c for c in out_cols if c in df.columns]]
    df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
