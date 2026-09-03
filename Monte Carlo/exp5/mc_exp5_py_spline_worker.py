#!/usr/bin/env python3
"""Run the spline-normalizer sensitivity variants for Experiment 5.

The worker reuses ``data.csv`` from Experiment 4 replications 1--240 and
changes one training dimension at a time. Architecture variants jointly vary
the embedding network, number of spline bins, and spline-parameter network.

Usage:
  python mc_exp5_py_spline_worker.py <task_id> <n_reps> <n_nodes> <n_cores>
    <source_output_dir> <variant_output_dir> <variant_name>
"""

import json
import math
import os
import shutil
import sys
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed


SCRIPT_DIR = Path(__file__).resolve().parent
EXP4_DIR = SCRIPT_DIR.parent / "exp4"
if not (EXP4_DIR / "mc_exp4_py_worker.py").exists():
    EXP4_DIR = SCRIPT_DIR
sys.path.insert(0, str(EXP4_DIR))

import mc_exp4_py_worker as base  # noqa: E402


RESULT_PREFIX = "medflow_spline"
MODEL_DIR_NAME = "models_spline"
PSE_PREFIX = "pse_spline"
INTV_PREFIX = "intv_spline"

DEFAULT_CONFIG = {
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

VARIANT_OVERRIDES = {
    "arch_up": {
        "dimension": "architecture",
        "direction": "up",
        "emb_net": [140, 120, 100, 80, 60, 40],
        "num_bins": 320,
        "spline_hidden_dims": (100, 80, 60, 40, 30, 20),
    },
    "arch_down": {
        "dimension": "architecture",
        "direction": "down",
        "emb_net": [80, 60, 40, 20],
        "num_bins": 192,
        "spline_hidden_dims": (40, 30, 20, 10),
    },
    "lr_up": {
        "dimension": "learning_rate",
        "direction": "up",
        "learning_rate": 3e-4,
    },
    "lr_down": {
        "dimension": "learning_rate",
        "direction": "down",
        "learning_rate": 3e-5,
    },
    "batch_up": {
        "dimension": "training_batch_size",
        "direction": "up",
        "trn_batch_size": 256,
    },
    "batch_down": {
        "dimension": "training_batch_size",
        "direction": "down",
        "trn_batch_size": 64,
    },
}

RESULT_COLUMNS = [
    "mf_ATE",
    "mf_DY",
    "mf_DM2Y",
    "mf_DM1Y",
    "mf_intv_OE",
    "mf_intv_IDE",
    "mf_intv_IIE",
]


def resolve_variant_config(variant_name):
    """Return the complete training configuration and named override."""
    if variant_name not in VARIANT_OVERRIDES:
        raise ValueError(
            f"Unknown variant '{variant_name}'. "
            f"Expected one of: {', '.join(sorted(VARIANT_OVERRIDES))}"
        )

    config = dict(DEFAULT_CONFIG)
    override = dict(VARIANT_OVERRIDES[variant_name])
    config.update({key: value for key, value in override.items() if key in DEFAULT_CONFIG})
    return config, override


def _validate_backend():
    info = base.report_runtime_backend("medflow-spline-exp5")
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
    required_arguments = ("norm_type", "num_bins", "bound", "spline_hidden_dims")
    missing_arguments = [name for name in required_arguments if name not in signature]
    if missing_arguments:
        raise RuntimeError(
            "The active cGNF.train signature is missing spline arguments: "
            + ", ".join(missing_arguments)
        )


def _result_path(output_dir, rep_id):
    result_dir = os.path.join(output_dir, "rep_results")
    os.makedirs(result_dir, exist_ok=True)
    return os.path.join(result_dir, f"{RESULT_PREFIX}_rep_{rep_id}.csv")


def _load_existing_success(output_dir, rep_id):
    result_path = _result_path(output_dir, rep_id)
    if not os.path.exists(result_path):
        return None

    try:
        result_df = pd.read_csv(result_path)
    except Exception:
        return None
    if result_df.empty:
        return None

    values = [
        pd.to_numeric(result_df.iloc[0][column], errors="coerce")
        for column in RESULT_COLUMNS
        if column in result_df.columns
    ]
    return result_df.iloc[0].to_dict() if any(pd.notna(value) for value in values) else None


def _ensure_rep_data(source_output_dir, variant_output_dir, rep_id):
    source_csv = os.path.join(source_output_dir, f"rep_{rep_id}", "data.csv")
    if not os.path.exists(source_csv):
        raise FileNotFoundError(f"{source_csv} not found")

    variant_rep_dir = os.path.join(variant_output_dir, f"rep_{rep_id}")
    os.makedirs(variant_rep_dir, exist_ok=True)
    variant_csv = os.path.join(variant_rep_dir, "data.csv")
    if not os.path.exists(variant_csv):
        shutil.copy2(source_csv, variant_csv)
    return variant_rep_dir


def _write_variant_metadata(
    variant_output_dir,
    source_output_dir,
    variant_name,
    variant_override,
    resolved_config,
    n_reps,
):
    metadata_path = os.path.join(variant_output_dir, "variant_config.json")
    payload = {
        "normalizer": "spline",
        "result_prefix": RESULT_PREFIX,
        "variant_name": variant_name,
        "source_output_dir": source_output_dir,
        "n_reps": int(n_reps),
        "default_config": DEFAULT_CONFIG,
        "override": variant_override,
        "resolved_config": resolved_config,
    }

    normalized_payload = json.loads(json.dumps(payload))
    if os.path.exists(metadata_path):
        with open(metadata_path, "r", encoding="utf-8") as handle:
            existing_payload = json.load(handle)
        protected_fields = (
            "normalizer",
            "result_prefix",
            "variant_name",
            "source_output_dir",
            "n_reps",
            "resolved_config",
        )
        mismatched_fields = [
            field
            for field in protected_fields
            if existing_payload.get(field) != normalized_payload.get(field)
        ]
        if mismatched_fields:
            raise RuntimeError(
                f"Existing metadata at {metadata_path} describes a different "
                f"run ({', '.join(mismatched_fields)} differ). Use a new output "
                "directory or remove the old run."
            )
        return

    temporary_path = f"{metadata_path}.{os.getpid()}.tmp"
    with open(temporary_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    os.replace(temporary_path, metadata_path)


def run_one_rep(rep_id, source_output_dir, variant_output_dir, variant_name, config):
    rep_id = int(rep_id)
    existing_result = _load_existing_success(variant_output_dir, rep_id)
    if existing_result is not None:
        print(f"[Rep {rep_id}] Existing successful spline result found. Skipping.")
        return existing_result

    _ensure_rep_data(source_output_dir, variant_output_dir, rep_id)
    result = base.run_one_rep(
        rep_id,
        variant_output_dir,
        train_kwargs=config,
        result_prefix=RESULT_PREFIX,
        model_dir_name=MODEL_DIR_NAME,
        pse_prefix=PSE_PREFIX,
        intv_prefix=INTV_PREFIX,
    )

    result["variant"] = variant_name
    result["normalizer"] = "spline"
    pd.DataFrame([result]).to_csv(_result_path(variant_output_dir, rep_id), index=False)
    return result


def main():
    if len(sys.argv) < 8:
        print(
            "Usage: python mc_exp5_py_spline_worker.py "
            "<task_id> <n_reps> <n_nodes> <n_cores> "
            "<source_output_dir> <variant_output_dir> <variant_name>"
        )
        return 1

    task_id = int(sys.argv[1])
    n_reps = int(sys.argv[2])
    n_nodes = int(sys.argv[3])
    n_cores = int(sys.argv[4])
    source_output_dir = sys.argv[5]
    variant_output_dir = sys.argv[6]
    variant_name = sys.argv[7]

    if n_reps != 240:
        raise ValueError(
            f"Experiment 5 must reuse exactly the first 240 replications; received {n_reps}."
        )
    if n_nodes < 1 or task_id < 1 or task_id > n_nodes:
        raise ValueError(
            f"task_id must be between 1 and n_nodes ({n_nodes}); received {task_id}."
        )

    _validate_backend()
    config, override = resolve_variant_config(variant_name)
    os.makedirs(variant_output_dir, exist_ok=True)
    _write_variant_metadata(
        variant_output_dir=variant_output_dir,
        source_output_dir=source_output_dir,
        variant_name=variant_name,
        variant_override=override,
        resolved_config=config,
        n_reps=n_reps,
    )

    reps_per_node = math.ceil(n_reps / n_nodes)
    start_rep = (task_id - 1) * reps_per_node + 1
    end_rep = min(task_id * reps_per_node, n_reps)
    replications = list(range(start_rep, end_rep + 1))

    print(f"Variant {variant_name}: {json.dumps(override, sort_keys=True)}")
    print(f"Resolved spline config: {json.dumps(config, sort_keys=True)}")
    print(
        f"Task {task_id} handling replications {start_rep} to {end_rep} "
        f"({len(replications)} total)"
    )
    print(f"Using {n_cores} parallel processes")

    Parallel(n_jobs=n_cores)(
        delayed(run_one_rep)(
            rep_id,
            source_output_dir,
            variant_output_dir,
            variant_name,
            config,
        )
        for rep_id in replications
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
