process RUN_PCA {
    label 'process_low'
    publishDir "${params.results_dir}/analysis", mode: 'copy'

    input:
    path(cohort_parquets)

    output:
    path("cross_cohort_pca.html"), emit: html

    script:
    """
    # Using the absolute path to ensure sklearn is found
    /home/nadeem/miniconda3/envs/variant-qc-pipeline/bin/python3 \
        ${projectDir}/bin/run_pca.py \
        --input-pq ${cohort_parquets} \
        --out-html cross_cohort_pca.html
    """
}
