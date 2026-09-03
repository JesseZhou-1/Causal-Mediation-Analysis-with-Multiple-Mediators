#!/usr/bin/env python3
"""
medflow spline worker for Experiment 4.

This wrapper uses the installed MedFlow package and requests the
spline-enabled cGNF backend through train_med(..., norm_type="spline").

Usage: python mc_exp4_py_spline_worker.py <task_id> <n_reps> <n_nodes> <n_cores> <output_dir>
"""

import math
import os
import sys
from pathlib import Path

from joblib import Parallel, delayed

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import mc_exp4_py_worker as base  # noqa: E402


TRAIN_KWARGS = {
    "trn_batch_size": 128,
    "val_batch_size": 2048,
    "learning_rate": 1e-4,
    "nb_epoch": 50000,
    "nb_estop": 200,
    "val_freq": 1,
    "emb_net": [100, 90, 80, 70, 60],
    "int_net": [60, 50, 40, 30, 20],
    "norm_type": "spline",
    "num_bins": 256,
    "bound": 7,
    "spline_hidden_dims": (60, 50, 40, 30, 20),
}


def _train_kwargs_with_env_overrides():
    train_kwargs = dict(TRAIN_KWARGS)
    env_map = {
        "MEDFLOW_NB_EPOCH": ("nb_epoch", int),
        "MEDFLOW_NB_ESTOP": ("nb_estop", int),
        "MEDFLOW_VAL_FREQ": ("val_freq", int),
        "MEDFLOW_TRN_BATCH_SIZE": ("trn_batch_size", int),
        "MEDFLOW_VAL_BATCH_SIZE": ("val_batch_size", int),
        "MEDFLOW_LEARNING_RATE": ("learning_rate", float),
    }
    for env_name, (key, cast) in env_map.items():
        value = os.environ.get(env_name)
        if value not in (None, ""):
            train_kwargs[key] = cast(value)
    return train_kwargs


def _validate_backend():
    info = base.report_runtime_backend("medflow-spline")
    if info.get("cgnf_error"):
        raise RuntimeError(f"Unable to import cGNF spline backend: {info['cgnf_error']}")
    if info.get("medflow_error"):
        raise RuntimeError(f"Unable to import MedFlow: {info['medflow_error']}")
    if not info.get("has_spline_normalizer", False):
        raise RuntimeError(
            "SplineNormalizer is not available in the active cGNF installation. "
            "Install the cGNF-spline backend before running this worker."
        )
    signature = info.get("cgnf_train_signature") or ""
    if "norm_type" not in signature:
        raise RuntimeError(
            "The active cGNF.train signature does not expose norm_type, so the "
            "spline backend is not active."
        )


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(
            "Usage: python mc_exp4_py_spline_worker.py "
            "<task_id> <n_reps> <n_nodes> <n_cores> <output_dir>"
        )
        sys.exit(1)

    _validate_backend()

    task_id = int(sys.argv[1])
    n_reps = int(sys.argv[2])
    n_nodes = int(sys.argv[3])
    n_cores = int(sys.argv[4])
    output_dir = sys.argv[5]

    reps_per_node = math.ceil(n_reps / n_nodes)
    start_rep = (task_id - 1) * reps_per_node + 1
    end_rep = min(task_id * reps_per_node, n_reps)
    my_reps = list(range(start_rep, end_rep + 1))

    print(
        f"Task {task_id} handling replications {start_rep} to {end_rep} "
        f"({len(my_reps)} total)"
    )
    print(f"Using {n_cores} parallel processes")
    train_kwargs = _train_kwargs_with_env_overrides()
    print(f"Spline train kwargs: {train_kwargs}")

    Parallel(n_jobs=n_cores)(
        delayed(base.run_one_rep)(
            rep_id,
            output_dir,
            train_kwargs=train_kwargs,
            result_prefix="medflow_spline",
            model_dir_name="models_spline",
            pse_prefix="pse_spline",
            intv_prefix="intv_spline",
        )
        for rep_id in my_reps
    )
