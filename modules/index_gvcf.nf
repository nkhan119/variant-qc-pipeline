// modules/index_gvcf.nf
// Index a gVCF with bcftools (tabix)

process INDEX_GVCF {
    label 'process_low'
    tag  "${cohort}/${sample}"

    input:
    tuple val(cohort), val(sample), path(gvcf)

    output:
    tuple val(cohort), val(sample), path(gvcf), path("${gvcf}.tbi")

    script:
    """
    bcftools index --tbi --threads ${task.cpus} ${gvcf}
    """
}
