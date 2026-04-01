// modules/variant_report.nf
// Generate per-cohort interactive Plotly HTML report + static PDF + parquet

process VARIANT_REPORT {
    label 'process_medium'
    tag  "${cohort}"

    publishDir "${params.results_dir}/${cohort}", mode: 'copy'

    input:
    tuple val(cohort), path(snv_tsv), path(cnv_tsv), path(metadata), path(multiqc_html)

    output:
    tuple val(cohort), path("${cohort}_variant_report.html"),     emit: html_report
    tuple val(cohort), path("${cohort}_variant_report.pdf"),      emit: pdf_report
    tuple val(cohort), path("${cohort}_snv_summary.parquet"),     emit: parquet
    tuple val(cohort), path("${cohort}_cnv_summary.parquet"),     emit: cnv_parquet
    tuple val(cohort), path("${cohort}_snv_table.tsv"),           emit: snv_table
    tuple val(cohort), path("${cohort}_cnv_table.tsv"),           emit: cnv_table

    script:
    """
    python3 ${projectDir}/bin/make_variant_report.py \
        --cohort        ${cohort} \
        --snv-tsv       ${snv_tsv} \
        --cnv-tsv       ${cnv_tsv} \
        --metadata      ${metadata} \
        --multiqc-html  ${multiqc_html} \
        --out-html      ${cohort}_variant_report.html \
        --out-pdf       ${cohort}_variant_report.pdf \
        --out-snv-pq    ${cohort}_snv_summary.parquet \
        --out-cnv-pq    ${cohort}_cnv_summary.parquet \
        --out-snv-tsv   ${cohort}_snv_table.tsv \
        --out-cnv-tsv   ${cohort}_cnv_table.tsv
    """
}
