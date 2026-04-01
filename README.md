# variant-qc-pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-≥23.04-brightgreen.svg)](https://www.nextflow.io/)
[![DSL2](https://img.shields.io/badge/DSL-2-blue.svg)](https://www.nextflow.io/docs/latest/dsl2.html)
[![Conda](https://img.shields.io/badge/conda-compatible-44A833.svg)](https://conda.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> A reproducible, performance-optimised Nextflow DSL2 pipeline for **SNV and CNV variant analysis** from gVCF cohorts, with MultiQC quality control, interactive **Plotly** HTML reports, and per-cohort PDF summaries.

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline DAG](#pipeline-dag)
3. [Features](#features)
4. [Requirements](#requirements)
5. [Installation](#installation)
6. [Input structure](#input-structure)
7. [Usage](#usage)
8. [Parameters](#parameters)
9. [Outputs](#outputs)
10. [Report figures](#report-figures)
11. [SLURM / Narval](#slurm--narval)
12. [Performance notes](#performance-notes)
13. [Author](#author)

---

## Overview

`variant-qc-pipeline` extends the [het-site-pipeline](https://github.com/nkhan119/het-site-pipeline2) with a full variant analysis layer:

- **MultiQC** integration over per-sample `bcftools stats` reports
- **SNV calling** — per-sample × per-chromosome scatter with hard-filter QC
- **CNV calling** — depth-based binning + ROH-proxy deletion calling
- **Annotation** — Ti/Tv class, AF approximation, het/hom flags per variant
- **Interactive Plotly HTML reports** — 15 figures per cohort (SNV QC, CNV analysis, genotype-phenotype)
- **Static PDF** exports via WeasyPrint
- **Parquet + TSV** summary tables for downstream R / Python analysis
- **Cross-cohort summary** — all cohorts merged into a single comparative report

---

## Pipeline DAG

```
GVCF input (per sample)
        │
        ▼
  INDEX_GVCF              (skip if .tbi already present)
        │
        ├──────────────────────────────────────────────┐
        ▼                                              ▼
  GVCF_STATS             CALL_SNVS × CALL_CNVS        │
  (bcftools stats)       (scatter: sample × chrom)    │
        │                        │                    │
        ▼                        ▼                    │
    MULTIQC              FILTER_SNVS                  │
  (per cohort)           (DP/GQ/QUAL filters)         │
        │                        │                    │
        │                        ▼                    │
        │               ANNOTATE_VARIANTS             │
        │               (Ti/Tv, AF, het flags)        │
        │                        │                    │
        └────────────────────────┤                    │
                                 ▼                    │
                           MERGE_COHORT ◄─────────────┘
                     (SNV + CNV + metadata)
                                 │
                                 ▼
                         VARIANT_REPORT
                   (Plotly HTML + PDF + parquet)
                                 │
                                 ▼
                         FINAL_SUMMARY
                  (cross-cohort HTML + parquet)
```

---

## Features

| Feature | Detail |
|---|---|
| Language | Nextflow DSL2 |
| SNV calling | `bcftools view` per chrom (scatter/gather) |
| SNV filters | DP ≥ 20, GQ ≥ 30, QUAL ≥ 30 (configurable) |
| CNV calling | Depth binning + bcftools ROH (configurable bin/min size) |
| QC | `bcftools stats` + MultiQC per cohort |
| Visualisation | **Plotly** interactive HTML (15 figures) + WeasyPrint PDF |
| Tables | Parquet + TSV per cohort + cross-cohort merged |
| Execution | Local · SLURM (Narval) · Docker (future) |
| Resumability | Full `-resume` support via Nextflow caching |
| Reports | HTML execution report · timeline · trace · DAG |

---

## Requirements

```
- Nextflow ≥ 23.04
- Conda / Mamba
- Java 11 or 17 (for Nextflow)
```

For SLURM on Narval:
```bash
module load StdEnv/2020
module load java/17
```

---

## Installation

```bash
# Clone the repository
git clone https://github.com/nkhan119/variant-qc-pipeline.git
cd variant-qc-pipeline

# Create and activate conda environment
conda env create -f environment.yaml
conda activate variant-qc-pipeline

# Verify Nextflow works
nextflow -version
```

---

## Input structure

Identical to het-site-pipeline — each cohort is a directory with gVCF files and a metadata TSV:

```
Cohort_A/
├── metadata.tsv              # Required columns: SampleID  (+ optional: Age, Ancestry, IQ, Sex)
├── sample1.gvcf.gz
├── sample1.gvcf.gz.tbi       # Auto-created if absent
├── sample2.gvcf.gz
└── ...
Cohort_B/
├── metadata.tsv
├── sampleX.gvcf.gz
└── ...
```

### metadata.tsv format

```
SampleID    Age    Ancestry    IQ    Sex
S001        32     EUR         108   F
S002        45     EAS         115   M
...
```

Only `SampleID` is mandatory. Optional columns (`Age`, `Ancestry`, `IQ`, `Sex`) enrich the phenotype-correlation plots.

---

## Usage

### Quick test (chr1 + chr22 only)

```bash
nextflow run main.nf -profile test
```

### Local full run

```bash
nextflow run main.nf -profile local

# Custom cohorts
nextflow run main.nf -profile local \
    --cohorts "Cohort_A,Cohort_B,Cohort_C"
```

### SLURM (Narval / Compute Canada)

```bash
sbatch submit_nextflow.sh

# Or with custom cohorts via environment variable
COHORTS="Cohort_A,Cohort_B" sbatch submit_nextflow.sh
```

### Resume after failure

```bash
nextflow run main.nf -profile slurm -resume
```

### Custom filters

```bash
nextflow run main.nf -profile slurm \
    --min_dp   30 \
    --min_gq   40 \
    --min_qual 50 \
    --cnv_bin_size 1000 \
    --cnv_min_size 5000
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `cohorts` | `Cohort_A,Cohort_B` | Comma-separated cohort directory names |
| `results_dir` | `results` | Root output directory |
| `min_dp` | `20` | Minimum FORMAT/DP for SNV filter |
| `min_gq` | `30` | Minimum FORMAT/GQ for SNV filter |
| `min_qual` | `30` | Minimum QUAL score for SNV filter |
| `min_af` | `0.05` | Minimum allele frequency (for annotation) |
| `cnv_bin_size` | `500` | CNV depth binning window (bp) |
| `cnv_min_size` | `1000` | Minimum CNV segment size to report (bp) |
| `chroms` | chr1–22,X,Y | Chromosomes for SNV scatter |
| `ref_genome` | `""` | Optional: path to reference FASTA |
| `max_memory` | `128.GB` | Resource cap |
| `max_cpus` | `32` | Resource cap |
| `max_time` | `24.h` | Resource cap |

All parameters can be overridden at runtime:

```bash
nextflow run main.nf -profile slurm --min_dp 30 --cohorts "Cohort_A,Cohort_C"
```

---

## Outputs

```
results/
├── Cohort_A/
│   ├── qc/
│   │   ├── stats/                      # Per-sample bcftools stats
│   │   ├── multiqc_Cohort_A.html       # MultiQC report
│   │   └── multiqc_Cohort_A_data/
│   ├── snv/
│   │   └── annotated/                  # Per-sample annotated SNV TSVs
│   ├── cnv/
│   │   └── raw/                        # Per-sample CNV BED files
│   ├── Cohort_A_snv_merged.tsv.gz      # All-sample SNV table
│   ├── Cohort_A_cnv_merged.tsv.gz      # All-sample CNV table
│   ├── Cohort_A_variant_report.html    # ★ Interactive Plotly report
│   ├── Cohort_A_variant_report.pdf     # ★ Static PDF report
│   ├── Cohort_A_snv_summary.parquet    # SNV summary (per-sample stats)
│   ├── Cohort_A_cnv_summary.parquet    # CNV summary (per-sample stats)
│   ├── Cohort_A_snv_table.tsv          # SNV summary TSV
│   └── Cohort_A_cnv_table.tsv          # CNV summary TSV
├── Cohort_B/
│   └── ...
├── summary/
│   ├── all_cohorts_snv.parquet         # Cross-cohort SNV parquet
│   ├── all_cohorts_cnv.parquet         # Cross-cohort CNV parquet
│   ├── cross_cohort_summary.html       # ★ Cross-cohort Plotly report
│   ├── cross_cohort_summary.pdf        # Cross-cohort PDF
│   ├── cross_cohort_snv_table.tsv
│   └── cross_cohort_cnv_table.tsv
└── reports/
    ├── execution_report.html            # Nextflow resource usage
    ├── timeline.html                    # Job timeline visualisation
    ├── trace.tsv                        # Per-task CPU/RAM/runtime
    └── pipeline_dag.html                # Pipeline DAG
```

---

## Report figures

### Per-cohort Plotly HTML report (15 figures)

**SNV QC**
1. SNV count per sample (bar + colour gradient)
2. Ti/Tv ratio per sample (bar + WGS reference line at 2.1)
3. SNV class breakdown — Ti / Tv / INDEL (stacked bar)
4. Allele frequency distribution (histogram + AF=0.5 marker)
5. Read depth (DP) distribution per sample (violin)
6. Genotype quality (GQ) distribution per sample (violin + GQ=30 filter line)
7. QUAL score distribution per sample (box)
8. Chromosomal SNV density (heatmap: samples × chromosomes)

**Genotype / Phenotype**
9. Heterozygous vs Hom-Alt counts coloured by ancestry (scatter)
10. SNV count vs Age (scatter + OLS regression)

**CNV Analysis**
11. CNV type breakdown per sample — GAIN / LOSS / DEL (stacked bar)
12. CNV size distribution log₁₀ by type (overlaid histogram)
13. Genome-wide log₂ copy-ratio per chromosome (strip plot)
14. CNV burden per sample (Mb affected) (bar)
15. CNV count by ancestry (box + points)

### Cross-cohort summary (8 figures)
Total SNV count, mean Ti/Tv, mean DP, het ratio, SNV violin, CNV count, CNV type breakdown, CNV burden per cohort.

---

## SLURM / Narval

The SLURM profile is pre-configured for Narval (Compute Canada):

| Setting | Value |
|---|---|
| Account | `def-fveyrier_cpu` |
| Queue | `cpu` |
| Queue size | 200 concurrent jobs |
| Poll interval | 30 seconds |
| Submit rate limit | 10 seconds |

Process resource labels:

| Label | CPUs | Memory | Time |
|---|---|---|---|
| `process_low` | 2 | 4 GB | 2 h |
| `process_medium` | 4 | 16 GB | 6 h |
| `process_high` | 8 | 32 GB | 12 h |

All jobs auto-retry (up to 2×) on memory/OOM exit codes (104, 134, 137, 139, 143, 247).

---

## Performance notes

- **SNV calling** is fully parallelised: each `(sample, chrom)` pair runs as an independent SLURM job → up to `n_samples × 24` concurrent tasks.
- **Nextflow `-resume`** skips any task whose inputs and process definition haven't changed — rerunning after a failure is instant.
- **CNV binning** uses pure Python + NumPy and is O(n_sites / bin_size) — fast even for WGS datasets.
- Parquet output (columnar, compressed) enables fast downstream analysis in R (`arrow`) or Python (`pandas` / `polars`).

---

## Citation / Reuse

If you use or adapt this pipeline, please credit the author.

##Author

Nadeem Khan, PhD Bioinformatician — INRS–Centre Armand-Frappier Santé-Biotechnologie, Laval, QC, Canada nkhan119@uottawa.ca @nkhan119

