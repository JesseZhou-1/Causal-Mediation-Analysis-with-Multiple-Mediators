#!/bin/bash
#SBATCH --job-name=mc_exp5_batch
#SBATCH --account=ssd
#SBATCH --time=24:00:00
#SBATCH --partition=ssd
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jessezhou1@rcc.uchicago.edu
#SBATCH --output=mc_exp5_batch_%j.out
#SBATCH --error=mc_exp5_batch_%j.err

if [ "$#" -lt 3 ]; then
    echo "Usage: sbatch run_mc_exp5_batch.sbatch <EXP4_OUTPUT_DIR> <EXP5_ROOT_DIR> <up|down>"
    exit 1
fi

SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$PWD}
if [ -f "$SUBMIT_DIR/medsim_mc/run_mc_exp5_group.sh" ]; then
    SCRIPT_DIR="$SUBMIT_DIR/medsim_mc"
elif [ -f "$SUBMIT_DIR/run_mc_exp5_group.sh" ]; then
    SCRIPT_DIR="$SUBMIT_DIR"
else
    echo "ERROR: Could not locate run_mc_exp5_group.sh from SLURM_SUBMIT_DIR=$SUBMIT_DIR"
    exit 1
fi

case "$3" in
    up|down)
        VARIANT="batch_$3"
        ;;
    *)
        echo "ERROR: third argument must be 'up' or 'down'"
        exit 1
        ;;
esac

bash "$SCRIPT_DIR/run_mc_exp5_group.sh" "$1" "$2" "$VARIANT"
