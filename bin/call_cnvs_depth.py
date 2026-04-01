#!/usr/bin/env python3
"""
call_cnvs_depth.py
------------------
Depth-based CNV detection from per-site DP values extracted from a gVCF.

Algorithm:
  1. Bin the genome into windows of --bin-size bp
  2. Compute median depth per bin
  3. Normalise by genome-wide median → log2 ratio
  4. Call gains (log2R > 0.58) and losses (log2R < -1.0) that span ≥ --min-size bp
  5. Merge adjacent bins of the same call type
  6. Append ROH segments (homozygous deletions) from bcftools roh output

Output: BED-like TSV with columns:
    CHROM  START  END  TYPE  LOG2R  N_BINS  SIZE  SAMPLE
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Chromosome sort order ──────────────────────────────────────
CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})


def chrom_sort_key(c):
    return CHROM_ORDER.get(c, 99)


def parse_args():
    p = argparse.ArgumentParser(description="Depth-based CNV caller")
    p.add_argument("--depth",    required=True,  help="Per-site depth TSV (CHROM POS DP)")
    p.add_argument("--sample",   required=True,  help="Sample ID")
    p.add_argument("--bin-size", type=int, default=500,  help="Bin size in bp")
    p.add_argument("--min-size", type=int, default=1000, help="Minimum CNV size in bp")
    p.add_argument("--roh",      default=None,   help="ROH TSV from bcftools roh (optional)")
    p.add_argument("--output",   required=True,  help="Output BED file")
    return p.parse_args()


def load_depth(depth_file: str) -> pd.DataFrame:
    df = pd.read_csv(
        depth_file, sep="\t", header=None,
        names=["CHROM", "POS", "DP"],
        dtype={"CHROM": str, "POS": int, "DP": float},
    )
    df = df[df["DP"] > 0].copy()
    return df


def bin_depth(df: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    """Aggregate depth into genomic bins."""
    df = df.copy()
    df["BIN"] = (df["POS"] // bin_size) * bin_size
    binned = (
        df.groupby(["CHROM", "BIN"])["DP"]
        .agg(["median", "count"])
        .reset_index()
        .rename(columns={"median": "MED_DP", "count": "N_SITES"})
    )
    binned = binned[binned["N_SITES"] >= 2].copy()
    return binned


def call_cnvs_from_bins(binned: pd.DataFrame, bin_size: int, min_size: int) -> pd.DataFrame:
    """Compute log2 ratios and segment into CNV calls."""
    genome_median = binned["MED_DP"].median()
    if genome_median <= 0:
        return pd.DataFrame()

    binned = binned.copy()
    binned["LOG2R"] = np.log2((binned["MED_DP"] + 0.5) / (genome_median + 0.5))

    # Classify bins
    binned["STATE"] = "NEUTRAL"
    binned.loc[binned["LOG2R"] >  0.58, "STATE"] = "GAIN"   # ~1.5× copy
    binned.loc[binned["LOG2R"] < -1.00, "STATE"] = "LOSS"   # ~0.5× copy (hetdel)
    binned.loc[binned["LOG2R"] < -3.00, "STATE"] = "DEL"    # homozygous deletion

    non_neutral = binned[binned["STATE"] != "NEUTRAL"].copy()
    if non_neutral.empty:
        return pd.DataFrame()

    # Sort by chrom order, then position
    non_neutral["_sort"] = non_neutral["CHROM"].map(chrom_sort_key)
    non_neutral = non_neutral.sort_values(["_sort", "BIN"])

    # Merge consecutive same-state bins on same chrom
    records = []
    for chrom, grp in non_neutral.groupby("CHROM", sort=False):
        grp = grp.sort_values("BIN")
        prev_state = None
        seg_start = seg_end = None
        seg_bins = []
        seg_log2rs = []

        for _, row in grp.iterrows():
            state = row["STATE"]
            pos   = row["BIN"]
            if state == prev_state and (pos - seg_end) <= bin_size * 2:
                seg_end = pos + bin_size
                seg_bins.append(1)
                seg_log2rs.append(row["LOG2R"])
            else:
                if prev_state and (seg_end - seg_start) >= min_size:
                    records.append({
                        "CHROM":  chrom,
                        "START":  seg_start,
                        "END":    seg_end,
                        "TYPE":   prev_state,
                        "LOG2R":  round(float(np.median(seg_log2rs)), 4),
                        "N_BINS": len(seg_bins),
                        "SIZE":   seg_end - seg_start,
                        "SOURCE": "DEPTH",
                    })
                prev_state = state
                seg_start  = pos
                seg_end    = pos + bin_size
                seg_bins   = [1]
                seg_log2rs = [row["LOG2R"]]

        # flush last segment
        if prev_state and seg_end is not None and (seg_end - seg_start) >= min_size:
            records.append({
                "CHROM":  chrom,
                "START":  seg_start,
                "END":    seg_end,
                "TYPE":   prev_state,
                "LOG2R":  round(float(np.median(seg_log2rs)), 4),
                "N_BINS": len(seg_bins),
                "SIZE":   seg_end - seg_start,
                "SOURCE": "DEPTH",
            })

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records)


def load_roh(roh_file: str, min_size: int) -> pd.DataFrame:
    """Parse bcftools roh output (RG lines) as candidate large deletions."""
    if not roh_file or not Path(roh_file).exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(
            roh_file, sep="\t", header=None,
            names=["CHROM", "START", "END", "TYPE", "QUAL"],
        )
        df["START"] = df["START"].astype(int)
        df["END"]   = df["END"].astype(int)
        df["SIZE"]  = df["END"] - df["START"]
        df = df[df["SIZE"] >= min_size].copy()
        df["LOG2R"]  = -2.0   # proxy for homozygous deletion
        df["N_BINS"] = 0
        df["SOURCE"] = "ROH"
        df["TYPE"]   = "DEL"
        return df[["CHROM", "START", "END", "TYPE", "LOG2R", "N_BINS", "SIZE", "SOURCE"]]
    except Exception as e:
        print(f"[WARN] Could not parse ROH file: {e}", file=sys.stderr)
        return pd.DataFrame()


def main():
    args = parse_args()

    print(f"[INFO] Loading depth file: {args.depth}", file=sys.stderr)
    df = load_depth(args.depth)

    if df.empty:
        print("[WARN] No depth data found. Writing empty output.", file=sys.stderr)
        header = "CHROM\tSTART\tEND\tTYPE\tLOG2R\tN_BINS\tSIZE\tSOURCE\tSAMPLE\n"
        Path(args.output).write_text(header)
        return

    print(f"[INFO] Binning {len(df):,} sites into {args.bin_size}bp windows", file=sys.stderr)
    binned  = bin_depth(df, args.bin_size)
    cnvs    = call_cnvs_from_bins(binned, args.bin_size, args.min_size)
    roh_df  = load_roh(args.roh, args.min_size)

    if not cnvs.empty and not roh_df.empty:
        combined = pd.concat([cnvs, roh_df], ignore_index=True)
    elif not cnvs.empty:
        combined = cnvs
    elif not roh_df.empty:
        combined = roh_df
    else:
        combined = pd.DataFrame(columns=["CHROM","START","END","TYPE","LOG2R","N_BINS","SIZE","SOURCE"])

    combined["SAMPLE"]  = args.sample
    combined["_sort"]   = combined["CHROM"].map(chrom_sort_key)
    combined = combined.sort_values(["_sort", "START"]).drop(columns=["_sort"])

    combined.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Called {len(combined)} CNV segments → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
