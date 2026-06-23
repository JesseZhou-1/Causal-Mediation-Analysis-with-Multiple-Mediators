#!/bin/bash
#SBATCH --job-name=mc_exp3
#SBATCH --account=ssd
#SBATCH --time=24:00:00
#SBATCH --partition=ssd
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jessezhou1@rcc.uchicago.edu
#SBATCH --output=mc_exp3_%j.out
#SBATCH --error=mc_exp3_%j.err

# ============================================================
# Monte Carlo Simulation for Experiment 3
# Count Treatment (Poisson), Binary M1, Continuous M2, Ordinal Y
# 1200 replications across 5 nodes x 48 cores
# ============================================================

echo "Job started at $(date)"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_NNODES"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

# Load modules
module load parallel
module load R/4.4.2+gcc-13.2.0

# Set R library path
export R_LIBS_USER=/scratch/midway3/jessezhou1/R_libs

# Working directory (where your R scripts are)
WORK_DIR=/project/wodtke/medsim_mc
cd $WORK_DIR

# Create output directory for this job
OUTPUT_DIR=${WORK_DIR}/results_exp3_${SLURM_JOB_ID}
mkdir -p $OUTPUT_DIR
echo "Output directory: $OUTPUT_DIR"

# Number of Monte Carlo replications
N_REPS=1200

# srun command for exclusive node access
srun="srun --exclusive -N1 -n1 --cpus-per-task=$SLURM_CPUS_PER_TASK"

# GNU parallel settings
parallel="parallel --delay 0.2 -j $SLURM_NNODES --joblog ${OUTPUT_DIR}/parallel_${SLURM_JOB_ID}.log --resume"

# Run Monte Carlo simulations
$parallel "$srun Rscript mc_exp3_worker.R {1} $N_REPS $SLURM_NNODES $SLURM_CPUS_PER_TASK $OUTPUT_DIR > ${OUTPUT_DIR}/task_{1}.log 2>&1" ::: $(seq 1 $SLURM_NNODES)

echo ""
echo "All tasks completed at $(date)"
echo "Aggregating results..."

# Aggregate results
Rscript aggregate_results_exp3.R $OUTPUT_DIR $N_REPS

echo "Job finished at $(date)"
