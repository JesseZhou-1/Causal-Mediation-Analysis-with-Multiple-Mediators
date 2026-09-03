#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: bash run_mc_exp5_spline_group.sh <EXP4_SPLINE_OUTPUT_DIR> <EXP5_SPLINE_ROOT_DIR> <variant_name>"
    exit 1
fi

EXP4_SPLINE_OUTPUT_DIR=$1
EXP5_SPLINE_ROOT_DIR=$2
VARIANT=$3
N_REPS=240
RESULT_PREFIX=medflow_spline

case "$VARIANT" in
    arch_up|arch_down)
        GROUP="arch"
        ;;
    lr_up|lr_down)
        GROUP="lr"
        ;;
    batch_up|batch_down)
        GROUP="batch"
        ;;
    *)
        echo "ERROR: variant_name must be one of: arch_up, arch_down, lr_up, lr_down, batch_up, batch_down"
        exit 1
        ;;
esac

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK_DIR=$SCRIPT_DIR

echo "Job started at $(date)"
echo "Experiment 5 spline group: $GROUP"
echo "Variant: $VARIANT"
echo "Source Experiment 4 spline directory: $EXP4_SPLINE_OUTPUT_DIR"
echo "Experiment 5 spline root directory: $EXP5_SPLINE_ROOT_DIR"
echo "Expected replications per variant: $N_REPS"
echo "Nodes: ${SLURM_NNODES:-NA}"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK:-NA}"

if [ ! -f "${EXP4_SPLINE_OUTPUT_DIR}/rep_1/data.csv" ]; then
    echo "ERROR: ${EXP4_SPLINE_OUTPUT_DIR}/rep_1/data.csv not found"
    exit 1
fi

if [ ! -f "${EXP4_SPLINE_OUTPUT_DIR}/rep_${N_REPS}/data.csv" ]; then
    echo "ERROR: ${EXP4_SPLINE_OUTPUT_DIR}/rep_${N_REPS}/data.csv not found"
    exit 1
fi

if [ ! -f "${EXP4_SPLINE_OUTPUT_DIR}/rep_results/${RESULT_PREFIX}_rep_1.csv" ] || \
   [ ! -f "${EXP4_SPLINE_OUTPUT_DIR}/rep_results/${RESULT_PREFIX}_rep_${N_REPS}.csv" ]; then
    echo "WARNING: baseline spline results for replications 1 and/or ${N_REPS} were not found."
    echo "The variants can run, but comparison against the spline baseline will be incomplete."
fi

if [ ! -f "${SCRIPT_DIR}/mc_exp5_py_spline_worker.py" ]; then
    echo "ERROR: ${SCRIPT_DIR}/mc_exp5_py_spline_worker.py not found"
    exit 1
fi

mkdir -p "$EXP5_SPLINE_ROOT_DIR"
cd "$WORK_DIR"

module load parallel
module load python/anaconda-2022.05

# Conda's activate/deactivate hooks are not nounset-safe.
set +u
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate /scratch/midway3/jessezhou1/myenv
set -u

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun_cmd="srun --exclusive -N1 -n1 --cpus-per-task=${SLURM_CPUS_PER_TASK}"
VARIANT_OUTPUT_DIR="${EXP5_SPLINE_ROOT_DIR}/${VARIANT}"
mkdir -p "$VARIANT_OUTPUT_DIR"

echo ""
echo "============================================================"
echo "Running spline variant: $VARIANT"
echo "Output directory: $VARIANT_OUTPUT_DIR"
echo "============================================================"

parallel_cmd="parallel --delay 0.2 -j ${SLURM_NNODES} --joblog ${VARIANT_OUTPUT_DIR}/parallel_py_${VARIANT}_${SLURM_JOB_ID}.log --resume"

$parallel_cmd "$srun_cmd python mc_exp5_py_spline_worker.py {1} $N_REPS ${SLURM_NNODES} ${SLURM_CPUS_PER_TASK} $EXP4_SPLINE_OUTPUT_DIR $VARIANT_OUTPUT_DIR $VARIANT > ${VARIANT_OUTPUT_DIR}/py_task_{1}.log 2>&1" ::: $(seq 1 ${SLURM_NNODES})

echo "Spline variant $VARIANT completed at $(date)"
echo "Aggregating spline results for $VARIANT..."

module load R/4.4.2+gcc-13.2.0
export R_LIBS_USER=/scratch/midway3/jessezhou1/R_libs
Rscript "$SCRIPT_DIR/aggregate_results_exp5_variant.R" \
    "$VARIANT_OUTPUT_DIR" "$N_REPS" "$VARIANT" "$RESULT_PREFIX"

echo ""
echo "Variant $VARIANT completed at $(date)"
echo "After all six variants finish, compare them with:"
echo "  Rscript $SCRIPT_DIR/compare_results_exp5.R $EXP4_SPLINE_OUTPUT_DIR $EXP5_SPLINE_ROOT_DIR $N_REPS $RESULT_PREFIX"
