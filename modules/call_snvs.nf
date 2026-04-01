// modules/call_snvs.nf
// Extract SNVs from a gVCF for a single chromosome using bcftools

process CALL_SNVS {
    label 'process_low'
    tag  "${cohort}/${sample}/${chrom}"

    input:
    tuple val(cohort), val(sample), path(gvcf), path(tbi), val(chrom)

    output:
    tuple val(cohort), val(sample), val(chrom), path("${sample}.${chrom}.raw_snv.vcf.gz"), emit: raw_vcf

    script:
    """
    bcftools view \
        --regions ${chrom} \
        --type snps \
        --min-ac 1 \
        --samples ${sample} \
        --output-type z \
        --threads ${task.cpus} \
        -o ${sample}.${chrom}.raw_snv.vcf.gz \
        ${gvcf}

    bcftools index --tbi ${sample}.${chrom}.raw_snv.vcf.gz
    """
}
