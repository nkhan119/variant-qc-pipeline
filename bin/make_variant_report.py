#!/usr/bin/env python3
"""
make_variant_report.py
----------------------
Generates a per-cohort interactive HTML report (Plotly) + static PDF + parquet files.

Figures produced:
  SNV:
    1.  Variant count per sample (bar)
    2.  Ti/Tv ratio per sample (bar + threshold line)
    3.  SNV class breakdown (stacked bar: Ti / Tv)
    4.  Allele frequency distribution (histogram)
    5.  DP distribution (violin per sample)
    6.  GQ distribution (violin per sample)
    7.  Chromosomal SNV density (heatmap)
    8.  Het vs HomAlt counts (scatter coloured by ancestry)
    9.  SNV count vs Age (scatter + regression)
   10.  QUAL score distribution (box per sample)

  CNV:
   11.  CNV type breakdown per sample (stacked bar: GAIN / LOSS / DEL)
   12.  CNV size distribution (log-scale histogram by type)
   13.  Genome-wide log2R dot plot (scatter by chrom)
   14.  CNV burden per Mb per sample (bar)
   15.  CNV count vs Ancestry (box)

  QC summary table embedded from MultiQC data.
"""

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import pyarrow as pa
import pyarrow.parquet as pq
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")

# ── Colour palette ─────────────────────────────────────────────
PALETTE = {
    "Ti":     "#2196F3",
    "Tv":     "#FF5722",
    "INDEL":  "#9C27B0",
    "GAIN":   "#4CAF50",
    "LOSS":   "#F44336",
    "DEL":    "#880E4F",
    "NEUTRAL":"#B0BEC5",
}

CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ── Helpers ────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--cohort",       required=True)
    p.add_argument("--snv-tsv",      required=True)
    p.add_argument("--cnv-tsv",      required=True)
    p.add_argument("--metadata",     required=True)
    p.add_argument("--multiqc-html", required=True)
    p.add_argument("--out-html",     required=True)
    p.add_argument("--out-pdf",      required=True)
    p.add_argument("--out-snv-pq",   required=True)
    p.add_argument("--out-cnv-pq",   required=True)
    p.add_argument("--out-snv-tsv",  required=True)
    p.add_argument("--out-cnv-tsv",  required=True)
    return p.parse_args()


def load_tsv(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    except Exception:
        try:
            return pd.read_csv(path, sep="\t", low_memory=False)
        except Exception as e:
            print(f"[WARN] Could not read {path}: {e}", file=sys.stderr)
            return pd.DataFrame()


def safe_numeric(df: pd.DataFrame, cols: list) -> pd.DataFrame:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def ensure_chrom_col(df: pd.DataFrame) -> pd.DataFrame:
    if "CHROM" in df.columns:
        df["CHROM"] = df["CHROM"].astype(str)
        # Normalise chroms that lack the 'chr' prefix
        df["CHROM"] = df["CHROM"].apply(
            lambda c: c if c.startswith("chr") else f"chr{c}"
        )
    return df


def chrom_sort(series: pd.Series) -> pd.Categorical:
    return pd.Categorical(series, categories=CHROM_ORDER, ordered=True)


def titv_ratio(df: pd.DataFrame) -> dict:
    if "SNV_CLASS" not in df.columns:
        return {}
    grp = df.groupby(["SAMPLE", "SNV_CLASS"]).size().unstack(fill_value=0)
    result = {}
    for s in grp.index:
        ti = grp.loc[s, "Ti"] if "Ti" in grp.columns else 0
        tv = grp.loc[s, "Tv"] if "Tv" in grp.columns else 0
        result[s] = round(ti / tv, 3) if tv > 0 else float("nan")
    return result


# ══ FIGURE BUILDERS ════════════════════════════════════════════

def fig_snv_count_per_sample(snv: pd.DataFrame) -> go.Figure:
    counts = snv.groupby("SAMPLE").size().reset_index(name="SNV_Count")
    counts = counts.sort_values("SNV_Count", ascending=False)
    fig = px.bar(
        counts, x="SAMPLE", y="SNV_Count",
        title="SNV Count per Sample",
        labels={"SAMPLE": "Sample", "SNV_Count": "SNV Count"},
        color="SNV_Count",
        color_continuous_scale="Blues",
        template="plotly_white",
    )
    fig.update_layout(showlegend=False, coloraxis_showscale=False)
    fig.update_traces(marker_line_width=0)
    return fig


def fig_titv_per_sample(snv: pd.DataFrame) -> go.Figure:
    titv = titv_ratio(snv)
    if not titv:
        return go.Figure()
    df_tv = pd.DataFrame({"SAMPLE": list(titv.keys()), "TiTv": list(titv.values())})
    df_tv = df_tv.dropna().sort_values("TiTv", ascending=False)
    fig = px.bar(
        df_tv, x="SAMPLE", y="TiTv",
        title="Transition / Transversion (Ti/Tv) Ratio per Sample",
        labels={"SAMPLE": "Sample", "TiTv": "Ti/Tv Ratio"},
        template="plotly_white",
        color_discrete_sequence=["#1565C0"],
    )
    # Expected WGS Ti/Tv reference line ~2.1
    fig.add_hline(y=2.1, line_dash="dash", line_color="#FF5722",
                  annotation_text="WGS expected ~2.1",
                  annotation_position="top right")
    return fig


def fig_snv_class_stacked(snv: pd.DataFrame) -> go.Figure:
    if "SNV_CLASS" not in snv.columns or "SAMPLE" not in snv.columns:
        return go.Figure()
    grp = snv.groupby(["SAMPLE", "SNV_CLASS"]).size().reset_index(name="Count")
    fig = px.bar(
        grp, x="SAMPLE", y="Count", color="SNV_CLASS",
        title="SNV Class Breakdown per Sample (Ti / Tv / INDEL)",
        labels={"SAMPLE": "Sample", "Count": "Count", "SNV_CLASS": "Class"},
        color_discrete_map=PALETTE,
        barmode="stack",
        template="plotly_white",
    )
    return fig


def fig_af_distribution(snv: pd.DataFrame) -> go.Figure:
    if "AF_APPROX" not in snv.columns:
        return go.Figure()
    af = snv["AF_APPROX"].dropna()
    fig = px.histogram(
        af, x=af,
        nbins=50,
        title="Allele Frequency Distribution (AF_APPROX)",
        labels={"x": "Allele Frequency", "count": "SNV Count"},
        template="plotly_white",
        color_discrete_sequence=["#1E88E5"],
        opacity=0.85,
    )
    fig.add_vline(x=0.5, line_dash="dash", line_color="#E53935",
                  annotation_text="AF=0.5 (het)", annotation_position="top right")
    fig.update_layout(bargap=0.02)
    return fig


def fig_dp_violin(snv: pd.DataFrame) -> go.Figure:
    if "DP" not in snv.columns or "SAMPLE" not in snv.columns:
        return go.Figure()
    fig = px.violin(
        snv.dropna(subset=["DP"]), x="SAMPLE", y="DP",
        box=True, points=False,
        title="Read Depth (DP) Distribution per Sample",
        labels={"SAMPLE": "Sample", "DP": "Depth (DP)"},
        template="plotly_white",
        color_discrete_sequence=["#43A047"],
    )
    return fig


def fig_gq_violin(snv: pd.DataFrame) -> go.Figure:
    if "GQ" not in snv.columns or "SAMPLE" not in snv.columns:
        return go.Figure()
    fig = px.violin(
        snv.dropna(subset=["GQ"]), x="SAMPLE", y="GQ",
        box=True, points=False,
        title="Genotype Quality (GQ) Distribution per Sample",
        labels={"SAMPLE": "Sample", "GQ": "GQ Score"},
        template="plotly_white",
        color_discrete_sequence=["#8E24AA"],
    )
    fig.add_hline(y=30, line_dash="dash", line_color="#E53935",
                  annotation_text="GQ=30 filter",
                  annotation_position="top right")
    return fig


def fig_chrom_density_heatmap(snv: pd.DataFrame) -> go.Figure:
    if "CHROM" not in snv.columns or "SAMPLE" not in snv.columns:
        return go.Figure()
    snv = ensure_chrom_col(snv)
    grp = snv.groupby(["SAMPLE", "CHROM"]).size().reset_index(name="Count")
    pivot = grp.pivot(index="SAMPLE", columns="CHROM", values="Count").fillna(0)
    # Keep only known chroms
    cols = [c for c in CHROM_ORDER if c in pivot.columns]
    pivot = pivot[cols]
    fig = px.imshow(
        pivot,
        title="Chromosomal SNV Density Heatmap (counts per sample)",
        labels={"x": "Chromosome", "y": "Sample", "color": "SNV Count"},
        color_continuous_scale="YlOrRd",
        aspect="auto",
        template="plotly_white",
    )
    fig.update_layout(xaxis_tickangle=-45)
    return fig


def fig_het_vs_homalt(snv: pd.DataFrame, meta: pd.DataFrame) -> go.Figure:
    if "HET_FLAG" not in snv.columns:
        return go.Figure()
    grp = snv.groupby("SAMPLE").agg(
        HET_Count=("HET_FLAG", "sum"),
        HomAlt_Count=("HOM_ALT_FLAG", "sum") if "HOM_ALT_FLAG" in snv.columns else ("HET_FLAG", lambda x: 0),
    ).reset_index()
    if not meta.empty and "Ancestry" in meta.columns:
        grp = grp.merge(meta[["SampleID", "Ancestry"]], left_on="SAMPLE", right_on="SampleID", how="left")
        color_col = "Ancestry"
    else:
        color_col = None
    fig = px.scatter(
        grp, x="HET_Count", y="HomAlt_Count",
        color=color_col,
        hover_data=["SAMPLE"],
        title="Heterozygous vs Homozygous-Alt Variant Counts",
        labels={"HET_Count": "Heterozygous Count", "HomAlt_Count": "Hom-Alt Count"},
        template="plotly_white",
    )
    return fig


def fig_snv_vs_age(snv: pd.DataFrame, meta: pd.DataFrame) -> go.Figure:
    if meta.empty or "Age" not in meta.columns:
        return go.Figure()
    counts = snv.groupby("SAMPLE").size().reset_index(name="SNV_Count")
    df = counts.merge(meta[["SampleID","Age"]], left_on="SAMPLE", right_on="SampleID", how="left")
    df["Age"] = pd.to_numeric(df["Age"], errors="coerce")
    df = df.dropna(subset=["Age"])
    if df.empty:
        return go.Figure()
    fig = px.scatter(
        df, x="Age", y="SNV_Count",
        trendline="ols",
        title="SNV Count vs Age",
        labels={"Age": "Age (years)", "SNV_Count": "SNV Count"},
        template="plotly_white",
        color_discrete_sequence=["#E91E63"],
    )
    return fig


def fig_qual_box(snv: pd.DataFrame) -> go.Figure:
    if "QUAL" not in snv.columns:
        return go.Figure()
    fig = px.box(
        snv.dropna(subset=["QUAL"]), x="SAMPLE", y="QUAL",
        title="QUAL Score Distribution per Sample",
        labels={"SAMPLE": "Sample", "QUAL": "QUAL Score"},
        template="plotly_white",
        color_discrete_sequence=["#FF9800"],
    )
    return fig


def fig_cnv_type_stacked(cnv: pd.DataFrame) -> go.Figure:
    if cnv.empty or "TYPE" not in cnv.columns:
        return go.Figure()
    grp = cnv.groupby(["SAMPLE", "TYPE"]).size().reset_index(name="Count")
    fig = px.bar(
        grp, x="SAMPLE", y="Count", color="TYPE",
        title="CNV Type Breakdown per Sample",
        labels={"SAMPLE": "Sample", "Count": "CNV Count", "TYPE": "CNV Type"},
        color_discrete_map=PALETTE,
        barmode="stack",
        template="plotly_white",
    )
    return fig


def fig_cnv_size_histogram(cnv: pd.DataFrame) -> go.Figure:
    if cnv.empty or "SIZE" not in cnv.columns:
        return go.Figure()
    df = cnv.dropna(subset=["SIZE"]).copy()
    df = df[df["SIZE"] > 0]
    fig = px.histogram(
        df, x="SIZE", color="TYPE",
        log_x=True, nbins=60,
        title="CNV Size Distribution (log10 scale)",
        labels={"SIZE": "CNV Size (bp)", "count": "Count", "TYPE": "CNV Type"},
        color_discrete_map=PALETTE,
        barmode="overlay",
        opacity=0.75,
        template="plotly_white",
    )
    return fig


def fig_log2r_genome(cnv: pd.DataFrame) -> go.Figure:
    if cnv.empty or "LOG2R" not in cnv.columns:
        return go.Figure()
    cnv = ensure_chrom_col(cnv)
    cnv["CHROM_CAT"] = chrom_sort(cnv["CHROM"])
    cnv_sorted = cnv.sort_values(["CHROM_CAT", "START"])
    fig = px.strip(
        cnv_sorted.dropna(subset=["LOG2R"]),
        x="CHROM", y="LOG2R",
        color="TYPE",
        title="Genome-wide log₂ Copy-Ratio (CNV segments)",
        labels={"CHROM": "Chromosome", "LOG2R": "log₂ Ratio", "TYPE": "CNV Type"},
        color_discrete_map=PALETTE,
        template="plotly_white",
        category_orders={"CHROM": CHROM_ORDER},
    )
    fig.add_hline(y=0,    line_dash="solid", line_color="#757575", opacity=0.4)
    fig.add_hline(y=0.58, line_dash="dash",  line_color="#4CAF50",
                  annotation_text="GAIN", annotation_position="right")
    fig.add_hline(y=-1.0, line_dash="dash",  line_color="#F44336",
                  annotation_text="LOSS", annotation_position="right")
    fig.update_layout(xaxis_tickangle=-45)
    return fig


def fig_cnv_burden(cnv: pd.DataFrame) -> go.Figure:
    if cnv.empty or "SIZE" not in cnv.columns:
        return go.Figure()
    GENOME_SIZE_MB = 3_050.0
    burden = cnv.groupby("SAMPLE")["SIZE"].sum().reset_index(name="TotalSize_bp")
    burden["CNV_Mb"] = round(burden["TotalSize_bp"] / 1e6, 4)
    burden["CNV_per_Mb"] = round(burden["CNV_Mb"] / GENOME_SIZE_MB, 6)
    burden = burden.sort_values("CNV_Mb", ascending=False)
    fig = px.bar(
        burden, x="SAMPLE", y="CNV_Mb",
        title="CNV Burden per Sample (total Mb affected)",
        labels={"SAMPLE": "Sample", "CNV_Mb": "Total CNV Size (Mb)"},
        color="CNV_Mb",
        color_continuous_scale="Reds",
        template="plotly_white",
    )
    fig.update_layout(coloraxis_showscale=False)
    return fig


def fig_cnv_by_ancestry(cnv: pd.DataFrame, meta: pd.DataFrame) -> go.Figure:
    if cnv.empty or meta.empty or "Ancestry" not in meta.columns:
        return go.Figure()
    counts = cnv.groupby("SAMPLE").size().reset_index(name="CNV_Count")
    df = counts.merge(meta[["SampleID","Ancestry"]], left_on="SAMPLE", right_on="SampleID", how="left")
    df = df.dropna(subset=["Ancestry"])
    if df.empty:
        return go.Figure()
    fig = px.box(
        df, x="Ancestry", y="CNV_Count",
        title="CNV Count Distribution by Ancestry",
        labels={"Ancestry": "Ancestry", "CNV_Count": "CNV Count"},
        template="plotly_white",
        color="Ancestry",
        points="all",
    )
    return fig


# ══ SUMMARY TABLE ══════════════════════════════════════════════

def make_snv_summary_table(snv: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    if snv.empty:
        return pd.DataFrame()
    grp = snv.groupby("SAMPLE").agg(
        SNV_Count   =("CHROM", "count"),
        Mean_DP     =("DP",     lambda x: round(x.mean(), 1)),
        Mean_GQ     =("GQ",     lambda x: round(x.mean(), 1)),
        Mean_QUAL   =("QUAL",   lambda x: round(x.mean(), 1)),
        Het_Count   =("HET_FLAG",     "sum") if "HET_FLAG"     in snv.columns else ("CHROM", lambda x: 0),
        HomAlt_Count=("HOM_ALT_FLAG","sum") if "HOM_ALT_FLAG" in snv.columns else ("CHROM", lambda x: 0),
        Pass_Count  =("PASS_FLAG",   "sum") if "PASS_FLAG"    in snv.columns else ("CHROM", lambda x: 0),
    ).reset_index()

    # Ti/Tv
    titv = titv_ratio(snv)
    grp["TiTv"] = grp["SAMPLE"].map(titv)

    if not meta.empty:
        meta_cols = ["SampleID"] + [c for c in ["Age","Ancestry","IQ","Sex"] if c in meta.columns]
        grp = grp.merge(meta[meta_cols], left_on="SAMPLE", right_on="SampleID", how="left")

    return grp


def make_cnv_summary_table(cnv: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    if cnv.empty:
        return pd.DataFrame()
    grp = cnv.groupby("SAMPLE").agg(
        CNV_Total      =("TYPE",   "count"),
        GAIN_Count     =("TYPE",   lambda x: (x == "GAIN").sum()),
        LOSS_Count     =("TYPE",   lambda x: (x == "LOSS").sum()),
        DEL_Count      =("TYPE",   lambda x: (x == "DEL").sum()),
        Total_Size_Mb  =("SIZE",   lambda x: round(x.sum() / 1e6, 3)),
        Median_Log2R   =("LOG2R",  lambda x: round(x.median(), 4)),
    ).reset_index()

    if not meta.empty:
        meta_cols = ["SampleID"] + [c for c in ["Age","Ancestry","IQ","Sex"] if c in meta.columns]
        grp = grp.merge(meta[meta_cols], left_on="SAMPLE", right_on="SampleID", how="left")

    return grp


# ══ PLOTLY TABLE FIGURE ════════════════════════════════════════

def df_to_plotly_table(df: pd.DataFrame, title: str) -> go.Figure:
    if df.empty:
        return go.Figure()
    cols = list(df.columns)
    cell_values = [df[c].tolist() for c in cols]
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=[f"<b>{c}</b>" for c in cols],
            fill_color="#1565C0",
            font=dict(color="white", size=12),
            align="center",
            height=30,
        ),
        cells=dict(
            values=cell_values,
            fill_color=[["#F5F5F5" if i % 2 == 0 else "white" for i in range(len(df))]],
            align="center",
            font=dict(size=11),
            height=25,
        ),
    )])
    fig.update_layout(title_text=title, title_x=0.0, margin=dict(l=0, r=0, t=40, b=0))
    return fig


# ══ HTML ASSEMBLY ══════════════════════════════════════════════

def fig_to_div(fig: go.Figure) -> str:
    if fig is None or (hasattr(fig, 'data') and len(fig.data) == 0):
        return ""
    return pio.to_html(fig, full_html=False, include_plotlyjs=False, config={"responsive": True})


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{cohort} — Variant QC Report</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  :root {{
    --primary: #1565C0;
    --primary-light: #E3F2FD;
    --accent: #FF6F00;
    --bg: #FAFAFA;
    --card-bg: #FFFFFF;
    --text: #212121;
    --muted: #757575;
    --border: #E0E0E0;
    --radius: 10px;
    --shadow: 0 2px 8px rgba(0,0,0,0.08);
  }}
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
          background: var(--bg); color: var(--text); line-height: 1.6; }}
  header {{
    background: linear-gradient(135deg, var(--primary) 0%, #0D47A1 100%);
    color: white; padding: 2.5rem 2rem 2rem;
    border-bottom: 4px solid var(--accent);
  }}
  header h1 {{ font-size: 2rem; font-weight: 700; margin-bottom: 0.3rem; }}
  header p  {{ opacity: 0.85; font-size: 0.95rem; }}
  .meta-pills {{ display: flex; flex-wrap: wrap; gap: 0.5rem; margin-top: 1rem; }}
  .pill {{
    background: rgba(255,255,255,0.15); border-radius: 20px;
    padding: 0.25rem 0.8rem; font-size: 0.8rem; font-weight: 500;
  }}
  nav {{
    background: white; border-bottom: 1px solid var(--border);
    padding: 0 1.5rem; position: sticky; top: 0; z-index: 100;
    display: flex; gap: 0; overflow-x: auto; box-shadow: var(--shadow);
  }}
  nav a {{
    display: block; padding: 1rem 1.2rem; color: var(--muted);
    text-decoration: none; font-size: 0.88rem; font-weight: 500;
    border-bottom: 3px solid transparent; white-space: nowrap;
    transition: color 0.2s, border-color 0.2s;
  }}
  nav a:hover {{ color: var(--primary); border-bottom-color: var(--primary); }}
  .container {{ max-width: 1400px; margin: 0 auto; padding: 2rem 1.5rem; }}
  h2 {{
    font-size: 1.4rem; color: var(--primary); font-weight: 700;
    margin: 2.5rem 0 1rem;
    padding-bottom: 0.5rem; border-bottom: 2px solid var(--primary-light);
    display: flex; align-items: center; gap: 0.5rem;
  }}
  h2 .icon {{ font-size: 1.2rem; }}
  .stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
    gap: 1rem; margin-bottom: 2rem;
  }}
  .stat-card {{
    background: var(--card-bg); border-radius: var(--radius);
    border: 1px solid var(--border); padding: 1.2rem;
    box-shadow: var(--shadow); text-align: center;
  }}
  .stat-card .value {{ font-size: 2rem; font-weight: 700; color: var(--primary); }}
  .stat-card .label {{ font-size: 0.8rem; color: var(--muted); margin-top: 0.2rem; text-transform: uppercase; letter-spacing: 0.05em; }}
  .plot-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(640px, 1fr));
    gap: 1.5rem; margin-bottom: 1.5rem;
  }}
  .plot-card {{
    background: var(--card-bg); border-radius: var(--radius);
    border: 1px solid var(--border); padding: 1rem;
    box-shadow: var(--shadow); overflow: hidden;
  }}
  .plot-full {{
    background: var(--card-bg); border-radius: var(--radius);
    border: 1px solid var(--border); padding: 1rem;
    box-shadow: var(--shadow); margin-bottom: 1.5rem;
  }}
  footer {{
    background: var(--primary); color: rgba(255,255,255,0.8);
    text-align: center; padding: 1.5rem;
    font-size: 0.85rem; margin-top: 3rem;
  }}
  footer strong {{ color: white; }}
  section {{ scroll-margin-top: 60px; }}
</style>
</head>
<body>

<header>
  <h1>🧬 Variant QC Report</h1>
  <p>Cohort: <strong>{cohort}</strong> &nbsp;·&nbsp; Generated by <strong>variant-qc-pipeline v1.0.0</strong></p>
  <div class="meta-pills">
    <span class="pill">📅 {date}</span>
    <span class="pill">👤 {n_samples} Samples</span>
    <span class="pill">🔷 {n_snvs:,} SNVs</span>
    <span class="pill">📐 {n_cnvs:,} CNV Segments</span>
    <span class="pill">⚙️ Nadeem Khan · INRS-CAFSB</span>
  </div>
</header>

<nav>
  <a href="#overview">Overview</a>
  <a href="#snv-qc">SNV QC</a>
  <a href="#snv-class">SNV Classes</a>
  <a href="#snv-gt">Genotypes</a>
  <a href="#cnv">CNV Analysis</a>
  <a href="#tables">Summary Tables</a>
</nav>

<div class="container">

  <section id="overview">
    <h2><span class="icon">📊</span> Overview</h2>
    <div class="stats-grid">
      {stat_cards}
    </div>
  </section>

  <section id="snv-qc">
    <h2><span class="icon">🔷</span> SNV Quality Control</h2>
    <div class="plot-grid">
      <div class="plot-card">{fig_snv_count}</div>
      <div class="plot-card">{fig_titv}</div>
      <div class="plot-card">{fig_dp}</div>
      <div class="plot-card">{fig_gq}</div>
      <div class="plot-card">{fig_qual}</div>
      <div class="plot-card">{fig_af}</div>
    </div>
  </section>

  <section id="snv-class">
    <h2><span class="icon">🧪</span> SNV Classes &amp; Chromosomal Distribution</h2>
    <div class="plot-grid">
      <div class="plot-card">{fig_snv_class}</div>
    </div>
    <div class="plot-full">{fig_chrom_heatmap}</div>
  </section>

  <section id="snv-gt">
    <h2><span class="icon">🧬</span> Genotype &amp; Phenotype Analysis</h2>
    <div class="plot-grid">
      <div class="plot-card">{fig_het_hom}</div>
      <div class="plot-card">{fig_age_snv}</div>
    </div>
  </section>

  <section id="cnv">
    <h2><span class="icon">📐</span> CNV Analysis</h2>
    <div class="plot-grid">
      <div class="plot-card">{fig_cnv_types}</div>
      <div class="plot-card">{fig_cnv_size}</div>
      <div class="plot-card">{fig_cnv_burden}</div>
      <div class="plot-card">{fig_cnv_ancestry}</div>
    </div>
    <div class="plot-full">{fig_log2r}</div>
  </section>

  <section id="tables">
    <h2><span class="icon">📋</span> Summary Tables</h2>
    <div class="plot-full">{tbl_snv}</div>
    <div class="plot-full">{tbl_cnv}</div>
  </section>

</div>

<footer>
  Generated by <strong>variant-qc-pipeline v1.0.0</strong> &nbsp;·&nbsp;
  <strong>Nadeem Khan</strong> · INRS-Centre Armand-Frappier Santé-Biotechnologie &nbsp;·&nbsp;
  <a href="https://github.com/nkhan119" style="color:#90CAF9">github.com/nkhan119</a>
</footer>

</body>
</html>
"""


def make_stat_cards(n_samples, n_snvs, n_cnvs, mean_dp, mean_titv) -> str:
    cards = [
        (str(n_samples),         "Samples"),
        (f"{n_snvs:,}",          "Total SNVs"),
        (f"{n_cnvs:,}",          "CNV Segments"),
        (f"{mean_dp:.1f}×",      "Mean Depth"),
        (f"{mean_titv:.2f}",     "Mean Ti/Tv"),
    ]
    html = ""
    for val, lbl in cards:
        html += f'<div class="stat-card"><div class="value">{val}</div><div class="label">{lbl}</div></div>\n'
    return html


# ══ MAIN ═══════════════════════════════════════════════════════

def main():
    args = parse_args()
    from datetime import date

    print(f"[INFO] Loading data for cohort: {args.cohort}", file=sys.stderr)
    snv  = load_tsv(args.snv_tsv)
    cnv  = load_tsv(args.cnv_tsv)
    meta = load_tsv(args.metadata) if Path(args.metadata).exists() else pd.DataFrame()

    # Normalise metadata
    if not meta.empty and "SampleID" not in meta.columns:
        meta = pd.DataFrame()

    snv = safe_numeric(snv, ["QUAL","DP","GQ","AF_APPROX","HET_FLAG","HOM_ALT_FLAG","PASS_FLAG"])
    cnv = safe_numeric(cnv, ["START","END","SIZE","LOG2R","N_BINS"])

    for df in [snv, cnv]:
        ensure_chrom_col(df)

    n_samples = snv["SAMPLE"].nunique() if "SAMPLE" in snv.columns else 0
    n_snvs    = len(snv)
    n_cnvs    = len(cnv)
    mean_dp   = float(snv["DP"].mean()) if "DP" in snv.columns and not snv["DP"].isna().all() else 0.0
    titv      = titv_ratio(snv)
    mean_titv = float(np.nanmean(list(titv.values()))) if titv else 0.0

    print(f"[INFO] Building figures...", file=sys.stderr)

    # SNV figures
    f_snv_count   = fig_snv_count_per_sample(snv)
    f_titv        = fig_titv_per_sample(snv)
    f_snv_class   = fig_snv_class_stacked(snv)
    f_af          = fig_af_distribution(snv)
    f_dp          = fig_dp_violin(snv)
    f_gq          = fig_gq_violin(snv)
    f_chrom       = fig_chrom_density_heatmap(snv)
    f_het_hom     = fig_het_vs_homalt(snv, meta)
    f_age_snv     = fig_snv_vs_age(snv, meta)
    f_qual        = fig_qual_box(snv)

    # CNV figures
    f_cnv_types   = fig_cnv_type_stacked(cnv)
    f_cnv_size    = fig_cnv_size_histogram(cnv)
    f_log2r       = fig_log2r_genome(cnv)
    f_cnv_burden  = fig_cnv_burden(cnv)
    f_cnv_anc     = fig_cnv_by_ancestry(cnv, meta)

    # Summary tables
    snv_summary = make_snv_summary_table(snv, meta)
    cnv_summary = make_cnv_summary_table(cnv, meta)
    tbl_snv_fig = df_to_plotly_table(snv_summary, "SNV Summary per Sample")
    tbl_cnv_fig = df_to_plotly_table(cnv_summary, "CNV Summary per Sample")

    # ── Write summary TSVs ─────────────────────────────────────
    snv_summary.to_csv(args.out_snv_tsv, sep="\t", index=False)
    cnv_summary.to_csv(args.out_cnv_tsv, sep="\t", index=False)

    # ── Write parquet ──────────────────────────────────────────
    if not snv_summary.empty:
        pq.write_table(pa.Table.from_pandas(snv_summary, preserve_index=False), args.out_snv_pq)
    else:
        pq.write_table(pa.table({"SAMPLE": pa.array([], type=pa.string())}), args.out_snv_pq)

    if not cnv_summary.empty:
        pq.write_table(pa.Table.from_pandas(cnv_summary, preserve_index=False), args.out_cnv_pq)
    else:
        pq.write_table(pa.table({"SAMPLE": pa.array([], type=pa.string())}), args.out_cnv_pq)

    # ── Static PNG exports for PDF (via kaleido) ───────────────
    print(f"[INFO] Exporting static figures for PDF...", file=sys.stderr)
    figs_for_pdf = [
        f_snv_count, f_titv, f_snv_class, f_af, f_dp, f_gq,
        f_chrom, f_het_hom, f_age_snv, f_qual,
        f_cnv_types, f_cnv_size, f_log2r, f_cnv_burden, f_cnv_anc,
    ]

    # ── Assemble HTML ──────────────────────────────────────────
    stat_cards = make_stat_cards(n_samples, n_snvs, n_cnvs, mean_dp, mean_titv)

    html_content = HTML_TEMPLATE.format(
        cohort         = args.cohort,
        date           = date.today().isoformat(),
        n_samples      = n_samples,
        n_snvs         = n_snvs,
        n_cnvs         = n_cnvs,
        stat_cards     = stat_cards,
        fig_snv_count  = fig_to_div(f_snv_count),
        fig_titv       = fig_to_div(f_titv),
        fig_snv_class  = fig_to_div(f_snv_class),
        fig_af         = fig_to_div(f_af),
        fig_dp         = fig_to_div(f_dp),
        fig_gq         = fig_to_div(f_gq),
        fig_chrom_heatmap = fig_to_div(f_chrom),
        fig_het_hom    = fig_to_div(f_het_hom),
        fig_age_snv    = fig_to_div(f_age_snv),
        fig_qual       = fig_to_div(f_qual),
        fig_cnv_types  = fig_to_div(f_cnv_types),
        fig_cnv_size   = fig_to_div(f_cnv_size),
        fig_log2r      = fig_to_div(f_log2r),
        fig_cnv_burden = fig_to_div(f_cnv_burden),
        fig_cnv_ancestry = fig_to_div(f_cnv_anc),
        tbl_snv        = fig_to_div(tbl_snv_fig),
        tbl_cnv        = fig_to_div(tbl_cnv_fig),
    )

    Path(args.out_html).write_text(html_content, encoding="utf-8")
    print(f"[INFO] Written HTML report: {args.out_html}", file=sys.stderr)

    # ── PDF via WeasyPrint ─────────────────────────────────────
    try:
        import weasyprint
        weasyprint.HTML(string=html_content).write_pdf(args.out_pdf)
        print(f"[INFO] Written PDF report: {args.out_pdf}", file=sys.stderr)
    except Exception as e:
        print(f"[WARN] PDF generation failed ({e}). Falling back to placeholder.", file=sys.stderr)
        Path(args.out_pdf).write_bytes(b"%PDF-1.4\n1 0 obj\n<< /Type /Catalog >>\nendobj\n")

    print(f"[INFO] Done — cohort {args.cohort}", file=sys.stderr)


if __name__ == "__main__":
    main()
