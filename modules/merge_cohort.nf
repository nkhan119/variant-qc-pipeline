// modules/merge_cohort.nf
// Merge all sample-level SNV and CNV tables into cohort-wide TSVs

process MERGE_COHORT {
    label 'process_medium'
    tag   "${cohort}"

    publishDir "${params.results_dir}/${cohort}", mode: 'copy'

    input:
    // We use underscores (_) for the sample ID lists because 
    // the python script doesn't actually need them; it just needs the files.
    tuple val(cohort), val(_), path(snv_tsvs), val(__), path(cnv_beds), path(metadata)

    output:
    tuple val(cohort), path("${cohort}_snv_merged.tsv.gz"), emit: cohort_snv
    tuple val(cohort), path("${cohort}_cnv_merged.tsv.gz"), emit: cohort_cnv

    script:
    """
    python3 ${projectDir}/bin/merge_cohort.py \
        --cohort   ${cohort} \
        --snv-dir  . \
        --cnv-dir  . \
        --metadata ${metadata} \
        --out-snv  ${cohort}_snv_merged.tsv.gz \
        --out-cnv  ${cohort}_cnv_merged.tsv.gz
    """
}
