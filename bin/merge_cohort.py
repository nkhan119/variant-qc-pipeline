#!/usr/bin/env python3
"""
merge_cohort.py
---------------
Merge all per-sample annotated SNV TSVs and CNV BED files into
cohort-wide compressed TSVs, joined with sample metadata.

Usage (called from Nextflow):
    python3 merge_cohort.py \\
        --cohort   Cohort_A \\
        --snv-dir  . \\
        --cnv-dir  . \\
        --metadata metadata.tsv \\
        --out-snv  Cohort_A_snv_merged.tsv.gz \\
        --out-cnv  Cohort_A_cnv_merged.tsv.gz
"""

import argparse
import glob
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--cohort",   required=True)
    p.add_argument("--snv-dir",  required=True)
    p.add_argument("--cnv-dir",  required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--out-snv",  required=True)
    p.add_argument("--out-cnv",  required=True)
    return p.parse_args()


def load_metadata(path: str) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str)
    if "SampleID" not in meta.columns:
        raise ValueError(f"metadata.tsv must contain a 'SampleID' column. Found: {list(meta.columns)}")
    return meta


def merge_snvs(snv_dir: str, meta: pd.DataFrame, cohort: str) -> pd.DataFrame:
    files = glob.glob(f"{snv_dir}/*.annotated.tsv.gz")
    if not files:
        print(f"[WARN] No SNV annotated files found in {snv_dir}", file=sys.stderr)
        return pd.DataFrame()

    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t", compression="gzip", dtype=str, low_memory=False)
            dfs.append(df)
        except Exception as e:
            print(f"[WARN] Could not read {f}: {e}", file=sys.stderr)

    if not dfs:
        return pd.DataFrame()

    merged = pd.concat(dfs, ignore_index=True)
    merged["COHORT"] = cohort

    # Join with metadata
    if "SAMPLE" in merged.columns and not meta.empty:
        merged = merged.merge(meta, left_on="SAMPLE", right_on="SampleID", how="left")

    # Cast numeric columns
    for col in ["QUAL", "DP", "GQ", "AF_APPROX"]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    print(f"[INFO] Merged SNVs: {len(merged):,} variants from {len(dfs)} samples", file=sys.stderr)
    return merged


def merge_cnvs(cnv_dir: str, meta: pd.DataFrame, cohort: str) -> pd.DataFrame:
    files = glob.glob(f"{cnv_dir}/*.cnv.bed")
    if not files:
        print(f"[WARN] No CNV BED files found in {cnv_dir}", file=sys.stderr)
        return pd.DataFrame()

    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
            dfs.append(df)
        except Exception as e:
            print(f"[WARN] Could not read {f}: {e}", file=sys.stderr)

    if not dfs:
        return pd.DataFrame()

    merged = pd.concat(dfs, ignore_index=True)
    merged["COHORT"] = cohort

    if "SAMPLE" in merged.columns and not meta.empty:
        merged = merged.merge(meta, left_on="SAMPLE", right_on="SampleID", how="left")

    for col in ["START", "END", "SIZE", "LOG2R", "N_BINS"]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    print(f"[INFO] Merged CNVs: {len(merged):,} segments from {len(dfs)} samples", file=sys.stderr)
    return merged


def main():
    args = parse_args()
    meta = load_metadata(args.metadata)

    snv_df = merge_snvs(args.snv_dir, meta, args.cohort)
    cnv_df = merge_cnvs(args.cnv_dir, meta, args.cohort)

    if snv_df.empty:
        snv_df = pd.DataFrame(columns=["SAMPLE","COHORT","CHROM","POS","REF","ALT","QUAL","DP","GQ"])

    if cnv_df.empty:
        cnv_df = pd.DataFrame(columns=["SAMPLE","COHORT","CHROM","START","END","TYPE","LOG2R","SIZE"])

    snv_df.to_csv(args.out_snv, sep="\t", index=False, compression="gzip")
    cnv_df.to_csv(args.out_cnv, sep="\t", index=False, compression="gzip")

    print(f"[INFO] Written: {args.out_snv}", file=sys.stderr)
    print(f"[INFO] Written: {args.out_cnv}", file=sys.stderr)


if __name__ == "__main__":
    main()
