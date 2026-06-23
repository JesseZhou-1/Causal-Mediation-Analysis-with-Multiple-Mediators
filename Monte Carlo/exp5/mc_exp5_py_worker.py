#!/usr/bin/env python3
"""
medflow-only worker for Experiment 5.
Reuses exp4 rep_{id}/data.csv for reps 1:240 and varies one hyperparameter
dimension at a time across six named variants.

Usage:
  python mc_exp5_py_worker.py <task_id> <n_reps> <n_nodes> <n_cores>
    <source_output_dir> <variant_output_dir> <variant_name>
"""

import json
import math
import os
import shutil
import sys
import traceback

import numpy as np
import pandas as pd

try:
    from medsim_mc.mc_exp4_py_worker import (
        CAT_LIST,
        _coerce_result_df,
        _extract_effects,
        _extract_intv_from_potential_outcomes,
        _extract_pse_from_potential_outcomes,
    )
except ImportError:
    from mc_exp4_py_worker import (
        CAT_LIST,
        _coerce_result_df,
        _extract_effects,
        _extract_intv_from_potential_outcomes,
        _extract_pse_from_potential_outcomes,
    )


DEFAULT_CONFIG = {
    "nb_epoch": 50000,
    "nb_estop": 50,
    "trn_batch_size": 128,
    "val_batch_size": 2048,
    "learning_rate": 1e-4,
    "emb_net": [100, 90, 80, 70, 60],
    "int_net": [60, 50, 40, 30, 20],
    "n_mce_samples": 5000,
}

VARIANT_OVERRIDES = {
    "arch_up": {
        "dimension": "architecture",
        "direction": "up",
        "emb_net": [140, 120, 100, 80, 60, 40],
        "int_net": [100, 80, 60, 40, 30, 20],
    },
    "arch_down": {
        "dimension": "architecture",
        "direction": "down",
        "emb_net": [80, 60, 40, 20],
        "int_net": [40, 30, 20, 10],
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

MF_RESULT_COLS = [
    "mf_ATE",
    "mf_DY",
    "mf_DM2Y",
    "mf_DM1Y",
    "mf_intv_OE",
    "mf_intv_IDE",
    "mf_intv_IIE",
]

PSE_PATTERNS = {
    "mf_ATE": r"TE\(|\bATE\b",
    "mf_DY": r"PSE.*D\s*[-=~]*>\s*Y|MNDE|direct",
    "mf_DM2Y": r"PSE.*M2\s*[-=~]*>\s*Y|via\s*M2",
    "mf_DM1Y": r"PSE.*M1\s*[-=~]*>\s*Y|via\s*M1",
}

INTV_PATTERNS = {
    "mf_intv_OE": r"OE\(",
    "mf_intv_IDE": r"IDE\(",
    "mf_intv_IIE": r"IIE\(",
}


def _resolve_variant_config(variant_name):
    if variant_name not in VARIANT_OVERRIDES:
        raise ValueError(
            f"Unknown variant '{variant_name}'. "
            f"Expected one of: {', '.join(sorted(VARIANT_OVERRIDES))}"
        )

    config = dict(DEFAULT_CONFIG)
    override = dict(VARIANT_OVERRIDES[variant_name])
    config.update({k: v for k, v in override.items() if k in DEFAULT_CONFIG})
    return config, override


def _rep_result_dir(output_dir):
    rep_result_dir = os.path.join(output_dir, "rep_results")
    os.makedirs(rep_result_dir, exist_ok=True)
    return rep_result_dir


def _result_path(output_dir, rep_id):
    return os.path.join(_rep_result_dir(output_dir), f"medflow_rep_{rep_id}.csv")


def _load_existing_success(output_dir, rep_id):
    out_path = _result_path(output_dir, rep_id)
    if not os.path.exists(out_path):
        return None

    try:
        df = pd.read_csv(out_path)
    except Exception:
        return None

    if df.empty:
        return None

    row = df.iloc[0].to_dict()
    numeric_vals = []
    for col in MF_RESULT_COLS:
        if col in df.columns:
            numeric_vals.append(pd.to_numeric(df.iloc[0][col], errors="coerce"))

    if any(pd.notna(v) for v in numeric_vals):
        return row

    return None


def _write_variant_metadata(variant_output_dir, source_output_dir, variant_name, variant_override, resolved_config, n_reps):
    meta_path = os.path.join(variant_output_dir, "variant_config.json")
    if os.path.exists(meta_path):
        return

    payload = {
        "variant_name": variant_name,
        "source_output_dir": source_output_dir,
        "n_reps": int(n_reps),
        "default_config": DEFAULT_CONFIG,
        "override": variant_override,
        "resolved_config": resolved_config,
    }

    try:
        with open(meta_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
    except Exception as exc:
        print(f"WARNING: failed to write {meta_path}: {exc}")


def _ensure_rep_data(source_output_dir, variant_output_dir, rep_id):
    source_rep_dir = os.path.join(source_output_dir, f"rep_{rep_id}")
    source_data_csv = os.path.join(source_rep_dir, "data.csv")
    if not os.path.exists(source_data_csv):
        raise FileNotFoundError(f"{source_data_csv} not found")

    variant_rep_dir = os.path.join(variant_output_dir, f"rep_{rep_id}")
    os.makedirs(variant_rep_dir, exist_ok=True)
    variant_data_csv = os.path.join(variant_rep_dir, "data.csv")
    if not os.path.exists(variant_data_csv):
        shutil.copy2(source_data_csv, variant_data_csv)

    return variant_rep_dir


def _save_result_row(output_dir, rep_id, row):
    out_path = _result_path(output_dir, rep_id)
    pd.DataFrame([row]).to_csv(out_path, index=False)
    print(f"[Rep {rep_id}] Results saved to {out_path}")


def run_one_rep(rep_id, source_output_dir, variant_output_dir, variant_name, variant_config):
    rep_id = int(rep_id)
    existing_row = _load_existing_success(variant_output_dir, rep_id)
    if existing_row is not None:
        print(f"[Rep {rep_id}] Existing successful result found. Skipping.")
        return existing_row

    results = {"rep_id": rep_id, "variant": variant_name}

    try:
        from medflow import sim_med, train_med

        rep_dir = _ensure_rep_data(source_output_dir, variant_output_dir, rep_id)
        data_csv = os.path.join(rep_dir, "data.csv")
        df = pd.read_csv(data_csv)
        print(f"[Rep {rep_id}] Data loaded: {df.shape[0]} rows, {df.shape[1]} columns")

        print(f"[Rep {rep_id}] Training medflow model for variant={variant_name}...")
        train_med(
            path=rep_dir + "/",
            dataset_name="data",
            treatment="D",
            mediator=["M1", "M2"],
            outcome="Y",
            confounder=["C"],
            cat_var=["M1", "Y"],
            nb_epoch=variant_config["nb_epoch"],
            nb_estop=variant_config["nb_estop"],
            trn_batch_size=variant_config["trn_batch_size"],
            val_batch_size=variant_config["val_batch_size"],
            learning_rate=variant_config["learning_rate"],
            seed=rep_id,
            emb_net=variant_config["emb_net"],
            int_net=variant_config["int_net"],
        )
        print(f"[Rep {rep_id}] Model training complete.")

        print(f"[Rep {rep_id}] Running PSE simulation...")
        pse_result = sim_med(
            path=rep_dir + "/",
            dataset_name="data",
            cat_list=CAT_LIST,
            intv_med=None,
            n_mce_samples=variant_config["n_mce_samples"],
            seed=rep_id,
            inv_datafile_name="pse",
        )
        print(f"[Rep {rep_id}] PSE simulation complete.")

        print(f"[Rep {rep_id}] Running interventional simulation...")
        intv_result = sim_med(
            path=rep_dir + "/",
            dataset_name="data",
            cat_list=CAT_LIST,
            intv_med=["M2=intv"],
            n_mce_samples=variant_config["n_mce_samples"],
            seed=rep_id + 100000,
            inv_datafile_name="intv",
        )
        intv_df = _coerce_result_df(intv_result, rep_dir, "intv")
        if intv_df is None or intv_df.empty:
            raise RuntimeError(f"[Rep {rep_id}] Missing interventional medflow summary output.")
        print(f"[Rep {rep_id}] Interventional simulation complete.")

        pse_df = _coerce_result_df(pse_result, rep_dir, "pse")
        if pse_df is None or pse_df.empty:
            raise RuntimeError(f"[Rep {rep_id}] Missing PSE medflow summary output.")
        if pse_df is not None:
            print(f"[Rep {rep_id}] PSE columns: {list(pse_df.columns)}")
        if intv_df is not None:
            print(f"[Rep {rep_id}] INTV columns: {list(intv_df.columns)}")

        pse_po = _extract_pse_from_potential_outcomes(pse_df, rep_id)
        intv_po = _extract_intv_from_potential_outcomes(intv_df, rep_id)

        if all(pd.isna(v) for v in pse_po.values()):
            print(f"[Rep {rep_id}] WARNING: PSE potential-outcome parsing failed; retrying with labeled regex extraction.")
            pse_po = _extract_effects(pse_df, PSE_PATTERNS, rep_id=rep_id, label="PSE")
        if all(pd.isna(v) for v in intv_po.values()):
            print(f"[Rep {rep_id}] WARNING: INTV potential-outcome parsing failed; retrying with labeled regex extraction.")
            intv_po = _extract_effects(intv_df, INTV_PATTERNS, rep_id=rep_id, label="INTV")

        if all(pd.isna(v) for v in pse_po.values()):
            raise RuntimeError(f"[Rep {rep_id}] Unable to extract PSE estimates from medflow output.")
        if all(pd.isna(v) for v in intv_po.values()):
            raise RuntimeError(f"[Rep {rep_id}] Unable to extract interventional estimates from medflow output.")

        results.update(pse_po)
        results.update(intv_po)

    except Exception as exc:
        err_dir = os.path.join(variant_output_dir, "error_logs")
        os.makedirs(err_dir, exist_ok=True)
        err_file = os.path.join(err_dir, f"medflow_rep_{rep_id}_error.log")
        with open(err_file, "w", encoding="utf-8") as handle:
            handle.write(f"rep_id: {rep_id}\n")
            handle.write(f"variant: {variant_name}\n")
            handle.write(f"error: {exc}\n\n")
            handle.write(traceback.format_exc())
        print(f"[Rep {rep_id}] ERROR: {exc}")
        print(f"[Rep {rep_id}] Error log written to {err_file}")

        for col in MF_RESULT_COLS:
            results[col] = np.nan

    _save_result_row(variant_output_dir, rep_id, results)
    return results


if __name__ == "__main__":
    if len(sys.argv) < 8:
        print(
            "Usage: python mc_exp5_py_worker.py <task_id> <n_reps> <n_nodes> <n_cores> "
            "<source_output_dir> <variant_output_dir> <variant_name>"
        )
        sys.exit(1)

    task_id = int(sys.argv[1])
    n_reps = int(sys.argv[2])
    n_nodes = int(sys.argv[3])
    n_cores = int(sys.argv[4])
    source_output_dir = sys.argv[5]
    variant_output_dir = sys.argv[6]
    variant_name = sys.argv[7]

    variant_config, variant_override = _resolve_variant_config(variant_name)
    os.makedirs(variant_output_dir, exist_ok=True)
    _write_variant_metadata(
        variant_output_dir=variant_output_dir,
        source_output_dir=source_output_dir,
        variant_name=variant_name,
        variant_override=variant_override,
        resolved_config=variant_config,
        n_reps=n_reps,
    )

    reps_per_node = math.ceil(n_reps / n_nodes)
    start_rep = (task_id - 1) * reps_per_node + 1
    end_rep = min(task_id * reps_per_node, n_reps)
    my_reps = list(range(start_rep, end_rep + 1))

    print(f"Variant {variant_name}: {json.dumps(variant_override, sort_keys=True)}")
    print(f"Task {task_id} handling replications {start_rep} to {end_rep} ({len(my_reps)} total)")
    print(f"Using {n_cores} parallel processes")

    from joblib import Parallel, delayed

    Parallel(n_jobs=n_cores)(
        delayed(run_one_rep)(
            rep_id,
            source_output_dir,
            variant_output_dir,
            variant_name,
            variant_config,
        )
        for rep_id in my_reps
    )
