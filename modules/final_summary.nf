// modules/final_summary.nf
// Merges all per-cohort parquets into a master cross-cohort dataset

process FINAL_SUMMARY {
    label 'process_medium'
    publishDir "${params.results_dir}/final_summary", mode: 'copy'

    input:
    path(snv_parquets)
    path(cnv_parquets)

    output:
    // These 'emit' labels are required so RUN_PCA can find the specific file
    path("all_cohorts_snv.parquet"),    emit: out_snv_pq
    path("all_cohorts_cnv.parquet"),    emit: out_cnv_pq
    path("cross_cohort_summary.html"),  emit: html
    path("cross_cohort_summary.pdf"),   emit: pdf
    path("cross_cohort_snv_table.tsv"), emit: snv_table
    path("cross_cohort_cnv_table.tsv"), emit: cnv_table

    script:
    """
    python3 ${projectDir}/bin/final_summary.py \
        --snv-parquets ${snv_parquets} \
        --cnv-parquets ${cnv_parquets} \
        --out-snv-pq    all_cohorts_snv.parquet \
        --out-cnv-pq    all_cohorts_cnv.parquet \
        --out-html      cross_cohort_summary.html \
        --out-pdf       cross_cohort_summary.pdf \
        --out-snv-tsv   cross_cohort_snv_table.tsv \
        --out-cnv-tsv   cross_cohort_cnv_table.tsv
    """
}
