// modules/filter_snvs.nf
// Apply hard filters to raw SNV calls

process FILTER_SNVS {
    label 'process_low'
    tag  "${cohort}/${sample}/${chrom}"

    input:
    tuple val(cohort), val(sample), val(chrom), path(raw_vcf)

    output:
    tuple val(cohort), val(sample), val(chrom), path("${sample}.${chrom}.filtered_snv.vcf.gz"), emit: filtered_vcf

    script:
    """
    # Hard filter: DP >= min_dp, GQ >= min_gq, QUAL >= min_qual
    bcftools filter \
        --include 'FORMAT/DP >= ${params.min_dp} && FORMAT/GQ >= ${params.min_gq} && QUAL >= ${params.min_qual}' \
        --output-type z \
        --threads ${task.cpus} \
        -o ${sample}.${chrom}.filtered_snv.vcf.gz \
        ${raw_vcf}

    bcftools index --tbi ${sample}.${chrom}.filtered_snv.vcf.gz
    """
}
