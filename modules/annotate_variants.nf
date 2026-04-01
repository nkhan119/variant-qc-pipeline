// modules/annotate_variants.nf
// Concatenate per-chrom VCFs and annotate SNVs to structured TSV

process ANNOTATE_VARIANTS {
    label 'process_medium'
    tag   "${cohort}/${sample}"

    publishDir "${params.results_dir}/${cohort}/snv/annotated", mode: 'copy', pattern: '*.annotated.tsv.gz'

    input:
    tuple val(cohort), val(sample), val(chroms), path(vcf_files)

    output:
    tuple val(cohort), val(sample), path("${sample}.annotated.tsv.gz"), emit: annotated

    script:
    // Ensure we have a space-separated string of files for the command line
    def vcf_list = vcf_files instanceof List ? vcf_files.join(' ') : vcf_files
    """
    # FIX: Generate indices for the staged VCF files on the fly
    # bcftools concat requires these to verify headers/positions
    for f in ${vcf_list}; do
        bcftools index -t \$f
    done

    # Concatenate all per-chrom filtered VCFs
    bcftools concat \
        --allow-overlaps \
        --remove-duplicates \
        --output-type z \
        --threads ${task.cpus} \
        -o ${sample}.merged_snv.vcf.gz \
        ${vcf_list}

    bcftools index --tbi ${sample}.merged_snv.vcf.gz

    # Extract rich annotation fields to TSV
    bcftools query \
        --samples ${sample} \
        --format '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t[%GT]\\t[%DP]\\t[%GQ]\\t[%AD]\\n' \
        ${sample}.merged_snv.vcf.gz \
    | awk -v OFS='\\t' 'BEGIN{print "CHROM","POS","ID","REF","ALT","QUAL","FILTER","GT","DP","GQ","AD"} {print}' \
    | python3 ${projectDir}/bin/annotate_snvs.py \
        --sample ${sample} \
        --cohort ${cohort} \
    | gzip -c > ${sample}.annotated.tsv.gz
    """
}
