#!/usr/bin/env nextflow

// ============================================================
//   variant-qc-pipeline  ·  Nextflow DSL2
//   Author : Nadeem Khan (nad119 · INRS-CAFSB)
//   Version: 1.1.1
// ============================================================

nextflow.enable.dsl = 2

// ── Parameters ───────────────────────────────────────────────
params.cohorts        = "Cohort_A,Cohort_B"
params.results_dir    = "results"
params.min_dp          = 20
params.min_gq          = 30
params.min_af          = 0.05
params.min_qual        = 30
params.cnv_min_size    = 1000
params.cnv_bin_size    = 500
params.ref_genome      = ""
params.chroms          = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10," +
                        "chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19," +
                        "chr20,chr21,chr22,chrX,chrY"

// ── Import modules ───────────────────────────────────────────
include { INDEX_GVCF          } from './modules/index_gvcf'
include { GVCF_STATS          } from './modules/gvcf_stats'
include { MULTIQC              } from './modules/multiqc'
include { CALL_SNVS           } from './modules/call_snvs'
include { FILTER_SNVS         } from './modules/filter_snvs'
include { CALL_CNVS           } from './modules/call_cnvs'
include { ANNOTATE_VARIANTS   } from './modules/annotate_variants'
include { MERGE_COHORT        } from './modules/merge_cohort'
include { VARIANT_REPORT      } from './modules/variant_report'
include { FINAL_SUMMARY       } from './modules/final_summary'
include { RUN_PCA             } from './modules/pca_analysis'

// ── Main workflow ─────────────────────────────────────────────
workflow {

    chrom_ch = Channel.from( params.chroms.split(",") )

    // Build per-sample channel: [cohort, sample, gvcf]
    gvcf_ch = Channel.from( params.cohorts.split(",") )
        .flatMap { cohort ->
            file("${cohort}/*.gvcf.gz").collect { gvcf ->
                def sample = gvcf.baseName.replaceAll(/\.gvcf$/, "")
                tuple( cohort, sample, gvcf )
            }
        }

    meta_ch = Channel.from( params.cohorts.split(",") )
        .map { cohort -> tuple( cohort, file("${cohort}/metadata.tsv") ) }

    // Step 1: Indexing
    gvcf_ch.branch {
        cohort, sample, gvcf ->
        needs_index : !file("${gvcf}.tbi").exists()
        has_index   : file("${gvcf}.tbi").exists()
    }.set { gvcf_branched }

    INDEX_GVCF( gvcf_branched.needs_index )
    
    indexed_ch = INDEX_GVCF.out.mix( 
        gvcf_branched.has_index.map { c, s, g -> tuple(c, s, g, file("${g}.tbi")) } 
    )

    // Step 2-3: QC & MultiQC
    GVCF_STATS( indexed_ch )
    MULTIQC( GVCF_STATS.out.stats_file.groupTuple( by: 0 ) )

    // Step 4-5: SNVs
    CALL_SNVS( indexed_ch.combine( chrom_ch ) )
    FILTER_SNVS( CALL_SNVS.out.raw_vcf )

    // Step 6: CNVs
    CALL_CNVS( indexed_ch )

    // Step 7: Annotation
    ANNOTATE_VARIANTS( FILTER_SNVS.out.filtered_vcf.groupTuple( by: [0, 1] ) )

    // Step 8: Merging
    snv_cohort = ANNOTATE_VARIANTS.out.annotated.groupTuple( by: 0 )
    cnv_cohort = CALL_CNVS.out.cnv_bed.groupTuple( by: 0 )
    
    merge_input = snv_cohort.join( cnv_cohort, by: 0 ).join( meta_ch, by: 0 )
    MERGE_COHORT( merge_input )

    // Step 9: Per-cohort variant reports
    report_input = MERGE_COHORT.out.cohort_snv
        .join( MERGE_COHORT.out.cohort_cnv, by: 0 )
        .join( meta_ch, by: 0 )
        .join( MULTIQC.out.report, by: 0 )

    VARIANT_REPORT( report_input )

    // ── Step 10: PCA Analysis (Exploratory) ──
    // FIX: Pulling from VARIANT_REPORT.out.parquet ensures we get actual Parquet files
    RUN_PCA( VARIANT_REPORT.out.parquet.map { it[1] }.collect() )

    // ── Step 11: Cross-cohort Final Summary (The Ultimate Step) ──
    FINAL_SUMMARY(
        VARIANT_REPORT.out.parquet.map { it[1] }.collect(),
        VARIANT_REPORT.out.cnv_parquet.map { it[1] }.collect()
    )
}
