#!/bin/bash
#SBATCH --job-name=mc_exp4_r
#SBATCH --account=ssd
#SBATCH --time=36:00:00
#SBATCH --partition=ssd
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jessezhou1@rcc.uchicago.edu
#SBATCH --output=mc_exp4_r_%j.out
#SBATCH --error=mc_exp4_r_%j.err

# ============================================================
# Experiment 4 — R Estimators (parametric methods)
# Nonlinear DGP: Count Treatment, Binary M1, Continuous M2, Ordinal Y
# Also saves data CSVs for medflow (run separately)
# ============================================================

echo "Job started at $(date)"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_NNODES"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"

# Working directory
WORK_DIR=/project/wodtke/medsim_mc
cd $WORK_DIR

# Create output directory for this job
OUTPUT_DIR=${WORK_DIR}/results_exp4_${SLURM_JOB_ID}
mkdir -p $OUTPUT_DIR
echo "Output directory: $OUTPUT_DIR"

# Number of Monte Carlo replications
N_REPS=600

# Load modules
module load parallel
module load R/4.4.2+gcc-13.2.0

# Set R library path
export R_LIBS_USER=/scratch/midway3/jessezhou1/R_libs

# Prevent ranger/OpenBLAS thread contention with mclapply workers
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# srun command for exclusive node access
srun="srun --exclusive -N1 -n1 --cpus-per-task=$SLURM_CPUS_PER_TASK"

# GNU parallel settings
parallel="parallel --delay 0.2 -j $SLURM_NNODES --joblog ${OUTPUT_DIR}/parallel_r_${SLURM_JOB_ID}.log --resume"

# Run R workers (generates data CSVs and runs R estimators)
$parallel "$srun Rscript mc_exp4_r_worker.R {1} $N_REPS $SLURM_NNODES $SLURM_CPUS_PER_TASK $OUTPUT_DIR > ${OUTPUT_DIR}/r_task_{1}.log 2>&1" ::: $(seq 1 $SLURM_NNODES)
PARALLEL_EXIT=$?

echo "GNU parallel exit code: ${PARALLEL_EXIT}"
if [[ ${PARALLEL_EXIT} -ne 0 ]]; then
  echo "Warning: one or more node tasks failed. Continuing to aggregation for partial recovery."
fi

JOBLOG_FILE="${OUTPUT_DIR}/parallel_r_${SLURM_JOB_ID}.log"
if [[ -f "${JOBLOG_FILE}" ]]; then
  echo "Failed task entries from GNU parallel joblog (Exitval != 0):"
  awk 'NR==1 || $7 != 0 {print}' "${JOBLOG_FILE}"
fi

MISSING_TASK_FILES=0
for task_id in $(seq 1 $SLURM_NNODES); do
  task_file="${OUTPUT_DIR}/results_task_${task_id}.rds"
  if [[ ! -f "${task_file}" ]]; then
    echo "Missing task result file: ${task_file}"
    MISSING_TASK_FILES=$((MISSING_TASK_FILES + 1))
  fi
done
echo "Missing task result files: ${MISSING_TASK_FILES}"

echo ""
echo "R estimators completed at $(date)"

# Aggregate R results
echo "Aggregating R results..."
Rscript aggregate_results_exp4.R $OUTPUT_DIR $N_REPS --phase=1
echo "Aggregation completed at $(date)"

echo "Job finished at $(date)"
