// modules/call_cnvs.nf
// Depth-based CNV detection using bcftools roh + custom Python binning

process CALL_CNVS {
    label 'process_medium'
    tag  "${cohort}/${sample}"

    publishDir "${params.results_dir}/${cohort}/cnv/raw", mode: 'copy', pattern: '*.cnv.bed'

    input:
    tuple val(cohort), val(sample), path(gvcf), path(tbi)

    output:
    tuple val(cohort), val(sample), path("${sample}.cnv.bed"), emit: cnv_bed

    script:
    """
    # Step 1: Extract per-site DP across all chroms
    bcftools query \
        --samples ${sample} \
        --format '%CHROM\t%POS\t[%DP]\n' \
        --regions ${params.chroms.replace(',', ',')} \
        ${gvcf} \
        | awk '\$3 != "." && \$3+0 > 0' \
        > ${sample}.depth.tsv

    # Step 2: Run ROH (runs of homozygosity) as proxy for large homozygous deletions
    bcftools roh \
        --AF-dflt 0.4 \
        --samples ${sample} \
        --output-type r \
        ${gvcf} \
        2>/dev/null \
        | grep -v "^#" \
        | awk -v OFS='\t' '{print \$3, \$4, \$5, "ROH", \$7}' \
        > ${sample}.roh.tsv || true

    # Step 3: Depth-based CNV calling using binning script
    python3 ${projectDir}/bin/call_cnvs_depth.py \
        --depth     ${sample}.depth.tsv \
        --sample    ${sample} \
        --bin-size  ${params.cnv_bin_size} \
        --min-size  ${params.cnv_min_size} \
        --roh       ${sample}.roh.tsv \
        --output    ${sample}.cnv.bed
    """
}
