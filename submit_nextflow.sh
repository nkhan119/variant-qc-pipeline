#!/bin/bash
#SBATCH --job-name=variant-qc
#SBATCH --account=def-fveyrier_cpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/nf_%j.out
#SBATCH --error=logs/nf_%j.err
#SBATCH --mail-type=FAIL

#Setup
set -e
mkdir -p logs

source ~/miniconda3/etc/profile.d/conda.sh
conda activate variant-qc-pipeline


nextflow run main.nf \
    -profile slurm \
    -resume \
    --cohorts "${COHORTS}" \
    --results_dir "results" \
    -with-report "results/reports/execution_report_${SLURM_JOB_ID}.html" \
    -with-trace "results/reports/trace_${SLURM_JOB_ID}.txt"

#Completion

EXIT_CODE=$?
echo "------------------------------------------------------------"
echo "  Finished at: $(date)"
echo "  Exit Code:   ${EXIT_CODE}"
echo "------------------------------------------------------------"

exit ${EXIT_CODE}
