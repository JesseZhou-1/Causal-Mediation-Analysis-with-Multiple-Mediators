#!/bin/bash
#SBATCH --job-name=mc_exp4_py_spline
#SBATCH --account=ssd
#SBATCH --time=36:00:00
#SBATCH --partition=ssd
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jessezhou1@rcc.uchicago.edu
#SBATCH --output=mc_exp4_py_spline_%j.out
#SBATCH --error=mc_exp4_py_spline_%j.err

# ============================================================
# Experiment 4 - MedFlow-RQS (Python) estimator
# Reads data CSVs from R job, trains cGNF-spline, estimates effects
#
# Usage: sbatch run_mc_exp4_py_spline.sbatch <R_OUTPUT_DIR>
#   e.g.: sbatch run_mc_exp4_py_spline.sbatch /project/wodtke/medsim_mc/results_exp4_spline
# ============================================================

# Get output directory from command line argument
if [ -z "$1" ]; then
    echo "ERROR: Must provide R job output directory as argument"
    echo "Usage: sbatch run_mc_exp4_py_spline.sbatch <R_OUTPUT_DIR>"
    exit 1
fi

OUTPUT_DIR=$1
N_REPS=600

echo "Job started at $(date)"
echo "Running on nodes: $SLURM_JOB_NODELIST"
echo "Number of nodes: $SLURM_NNODES"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Output directory: $OUTPUT_DIR"

# Verify R results exist
if [ ! -f "${OUTPUT_DIR}/r_results_exp4.rds" ]; then
    echo "ERROR: R results not found at ${OUTPUT_DIR}/r_results_exp4.rds"
    echo "Run the R job first (run_mc_exp4_r.sbatch)"
    exit 1
fi

# Locate the experiment scripts from the submission directory.
SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$PWD}
if [ -f "$SUBMIT_DIR/mc_exp4_py_spline_worker.py" ]; then
    WORK_DIR="$SUBMIT_DIR"
elif [ -f "$SUBMIT_DIR/Monte Carlo/exp4/mc_exp4_py_spline_worker.py" ]; then
    WORK_DIR="$SUBMIT_DIR/Monte Carlo/exp4"
elif [ -f "$SUBMIT_DIR/medsim_mc/mc_exp4_py_spline_worker.py" ]; then
    WORK_DIR="$SUBMIT_DIR/medsim_mc"
else
    echo "ERROR: Could not locate mc_exp4_py_spline_worker.py from $SUBMIT_DIR"
    exit 1
fi
cd "$WORK_DIR"

# Load Python environment
module load parallel

module load python/anaconda-2022.05

source activate

cd /scratch/midway3/jessezhou1/

conda activate /scratch/midway3/jessezhou1/myenv

cd "$WORK_DIR"

# Prevent PyTorch from using excessive threads
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# srun command for exclusive node access
srun="srun --exclusive -N1 -n1 --cpus-per-task=$SLURM_CPUS_PER_TASK"

# GNU parallel settings
parallel="parallel --delay 0.2 -j $SLURM_NNODES --joblog ${OUTPUT_DIR}/parallel_py_${SLURM_JOB_ID}.log --resume"

# Run medflow spline workers (same pattern as R job: one task per node, 48 cores each)
$parallel "$srun python mc_exp4_py_spline_worker.py {1} $N_REPS $SLURM_NNODES $SLURM_CPUS_PER_TASK $OUTPUT_DIR > ${OUTPUT_DIR}/py_spline_task_{1}.log 2>&1" ::: $(seq 1 $SLURM_NNODES)

echo ""
echo "medflow spline completed at $(date)"

# Final aggregation (R + medflow spline)
echo "Aggregating final results (R + medflow spline)..."
module load R/4.4.2+gcc-13.2.0
export R_LIBS_USER=/scratch/midway3/jessezhou1/R_libs
Rscript "$WORK_DIR/aggregate_results_exp4_spline.R" "$OUTPUT_DIR" "$N_REPS" --phase=2
echo "Aggregation completed at $(date)"

echo "Job finished at $(date)"
