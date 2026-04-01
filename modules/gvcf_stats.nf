// modules/gvcf_stats.nf
// Generate per-sample bcftools stats for MultiQC consumption

process GVCF_STATS {
    label 'process_low'
    tag  "${cohort}/${sample}"

    publishDir "${params.results_dir}/${cohort}/qc/stats", mode: 'copy', pattern: '*.stats'

    input:
    tuple val(cohort), val(sample), path(gvcf), path(tbi)

    output:
    tuple val(cohort), val(sample), path("${sample}.stats"), emit: stats_file

    script:
    """
    bcftools stats \
        --threads ${task.cpus} \
        --samples ${sample} \
        ${gvcf} \
        > ${sample}.stats
    """
}
