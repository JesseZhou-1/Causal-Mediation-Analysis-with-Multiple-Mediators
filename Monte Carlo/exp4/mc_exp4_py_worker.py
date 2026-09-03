#!/usr/bin/env python3
"""
medflow worker for Experiment 4.
Trains cGNF model, runs PSE and interventional simulations.

Usage: python mc_exp4_py_worker.py <task_id> <n_reps> <n_nodes> <n_cores> <output_dir>

Each task handles a subset of replications (same interface as R worker).
Expects data CSV at: {output_dir}/rep_{rep_id}/data.csv (written by R worker)
Saves results to:    {output_dir}/rep_results/medflow_rep_{rep_id}.csv
"""

import sys
import os
import traceback
import inspect
import pandas as pd
import numpy as np
import re

CONTRAST_DSTAR = 0.0
CONTRAST_D = 5.0
CAT_LIST = [int(CONTRAST_DSTAR), int(CONTRAST_D)]


def _rep_result_dir(output_dir):
    rep_result_dir = os.path.join(output_dir, "rep_results")
    os.makedirs(rep_result_dir, exist_ok=True)
    return rep_result_dir


def _runtime_backend_info():
    info = {
        "medflow_path": None,
        "medflow_error": None,
        "cgnf_path": None,
        "cgnf_train_signature": None,
        "cgnf_error": None,
        "has_spline_normalizer": False,
        "spline_error": None,
    }

    try:
        import medflow
        info["medflow_path"] = getattr(medflow, "__file__", None)
    except Exception as exc:
        info["medflow_error"] = f"{type(exc).__name__}: {exc}"

    try:
        import cGNF
        from cGNF import train as cgnf_train

        info["cgnf_path"] = getattr(cGNF, "__file__", None)
        info["cgnf_train_signature"] = str(inspect.signature(cgnf_train))
    except Exception as exc:
        info["cgnf_error"] = f"{type(exc).__name__}: {exc}"

    try:
        from cGNF.GNF_Modules.Normalizers import SplineNormalizer  # noqa: F401
        info["has_spline_normalizer"] = True
    except Exception as exc:
        info["has_spline_normalizer"] = False
        info["spline_error"] = f"{type(exc).__name__}: {exc}"

    return info


def report_runtime_backend(label="medflow"):
    info = _runtime_backend_info()
    print(f"[{label}] medflow path: {info['medflow_path']}")
    print(f"[{label}] medflow import error: {info['medflow_error']}")
    print(f"[{label}] cGNF path: {info['cgnf_path']}")
    print(f"[{label}] cGNF train signature: {info['cgnf_train_signature']}")
    print(f"[{label}] cGNF import error: {info['cgnf_error']}")
    print(f"[{label}] SplineNormalizer available: {info['has_spline_normalizer']}")
    print(f"[{label}] SplineNormalizer import error: {info['spline_error']}")
    return info


def _coerce_result_df(result_obj, rep_dir, inv_prefix):
    """Normalize medflow summary output to a pandas DataFrame."""
    if isinstance(result_obj, pd.DataFrame):
        return result_obj

    if isinstance(result_obj, str) and os.path.exists(result_obj):
        try:
            return pd.read_csv(result_obj)
        except Exception:
            pass

    if isinstance(result_obj, dict):
        try:
            return pd.DataFrame(result_obj)
        except Exception:
            pass

    if isinstance(result_obj, (list, tuple)):
        try:
            return pd.DataFrame(result_obj)
        except Exception:
            pass

    prefixes = [inv_prefix] if isinstance(inv_prefix, str) else list(inv_prefix)
    for pref in prefixes:
        candidate = os.path.join(rep_dir, f"{pref}_results.csv")
        if os.path.exists(candidate):
            try:
                return pd.read_csv(candidate)
            except Exception:
                return None
    return None


def _infer_value_col(df):
    if "Value" in df.columns:
        return "Value"
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if numeric_cols:
        return numeric_cols[-1]
    return None


def _extract_effects(df, effect_patterns, rep_id, label):
    """Extract scalar effects by regex over all non-value columns."""
    out = {k: np.nan for k in effect_patterns.keys()}

    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        print(f"[Rep {rep_id}] WARNING: {label} result is empty/None.")
        return out

    value_col = _infer_value_col(df)
    if value_col is None:
        print(f"[Rep {rep_id}] WARNING: {label} has no numeric/Value column. Columns={list(df.columns)}")
        return out

    text_cols = [c for c in df.columns if c != value_col]
    if text_cols:
        row_text = df[text_cols].astype(str).agg(" ".join, axis=1)
    else:
        row_text = pd.Series([""] * len(df))

    value_series = pd.to_numeric(df[value_col], errors="coerce")

    for key, pattern in effect_patterns.items():
        mask = row_text.str.contains(pattern, regex=True, case=False, na=False)
        if mask.any():
            out[key] = value_series.loc[mask].iloc[0]

    return out


def _get_label_col(df):
    for col in df.columns:
        lc = str(col).lower()
        if "potential" in lc and "outcome" in lc:
            return col
    for col in df.columns:
        if not pd.api.types.is_numeric_dtype(df[col]):
            return col
    return None


def _safe_float(x):
    try:
        return float(str(x).strip())
    except Exception:
        return None


def _parse_potential_rows(df):
    if df is None or not isinstance(df, pd.DataFrame) or df.empty:
        return [], None

    value_col = _infer_value_col(df)
    label_col = _get_label_col(df)
    if value_col is None or label_col is None:
        return [], value_col

    rows = []
    for _, row in df.iterrows():
        label = str(row[label_col])
        value = pd.to_numeric(row[value_col], errors="coerce")

        main_d = None
        m = re.search(r"Y\(\s*D\s*=\s*([^\s,\)]+)", label, flags=re.IGNORECASE)
        if m:
            main_d = m.group(1)

        m1_d = None
        m = re.search(r"M1\(\s*D\s*=\s*([^\s,\)]+)\)", label, flags=re.IGNORECASE)
        if m:
            m1_d = m.group(1)

        m2_d = None
        m = re.search(r"M2\(\s*D\s*=\s*([^\s,\)]+)\)", label, flags=re.IGNORECASE)
        if m:
            m2_d = m.group(1)

        m2_intv = bool(re.search(r"M2\s*=\s*intv", label, flags=re.IGNORECASE))

        rows.append({
            "label": label,
            "value": value,
            "main_d": main_d,
            "m1_d": m1_d,
            "m2_d": m2_d,
            "m2_intv": m2_intv,
        })
    return rows, value_col


def _find_value(rows, predicate):
    for r in rows:
        if predicate(r):
            return r["value"]
    return np.nan


def _same_level(parsed_level, target_level):
    parsed_float = _safe_float(parsed_level)
    if parsed_float is None:
        return False
    return parsed_float == float(target_level)


def _extract_pse_from_potential_outcomes(df, rep_id):
    out = {"mf_ATE": np.nan, "mf_DY": np.nan, "mf_DM2Y": np.nan, "mf_DM1Y": np.nan}
    rows, _ = _parse_potential_rows(df)
    if not rows:
        return out

    y_dstar = _find_value(
        rows,
        lambda r: _same_level(r["main_d"], CONTRAST_DSTAR) and r["m1_d"] is None and r["m2_d"] is None and not r["m2_intv"],
    )
    y_d = _find_value(
        rows,
        lambda r: _same_level(r["main_d"], CONTRAST_D) and r["m1_d"] is None and r["m2_d"] is None and not r["m2_intv"],
    )
    y_d_m1dstar = _find_value(
        rows,
        lambda r: _same_level(r["main_d"], CONTRAST_D) and _same_level(r["m1_d"], CONTRAST_DSTAR) and r["m2_d"] is None and not r["m2_intv"],
    )
    y_d_m1dstar_m2dstar = _find_value(
        rows,
        lambda r: _same_level(r["main_d"], CONTRAST_D) and _same_level(r["m1_d"], CONTRAST_DSTAR) and _same_level(r["m2_d"], CONTRAST_DSTAR) and not r["m2_intv"],
    )

    # Primary decomposition (aligns with D=d pathway formulas used in R truth):
    # ATE  = E[Y(d)] - E[Y(d*)]
    # DY   = E[Y(d, M1(d*), M2(d*))] - E[Y(d*)]
    # DM2Y = E[Y(d, M1(d*), M2(d))] - E[Y(d, M1(d*), M2(d*))]
    # DM1Y = E[Y(d)] - E[Y(d, M1(d*), M2(d))]
    if not any(pd.isna(v) for v in [y_dstar, y_d, y_d_m1dstar, y_d_m1dstar_m2dstar]):
        out["mf_ATE"] = y_d - y_dstar
        out["mf_DY"] = y_d_m1dstar_m2dstar - y_dstar
        out["mf_DM2Y"] = y_d_m1dstar - y_d_m1dstar_m2dstar
        out["mf_DM1Y"] = y_d - y_d_m1dstar
    return out


def _extract_intv_from_potential_outcomes(df, rep_id):
    out = {"mf_intv_OE": np.nan, "mf_intv_IDE": np.nan, "mf_intv_IIE": np.nan}
    rows, _ = _parse_potential_rows(df)
    if not rows:
        return out

    y_dstar = _find_value(
        rows,
        lambda r: _same_level(r["main_d"], CONTRAST_DSTAR) and r["m1_d"] is None and r["m2_d"] is None and not r["m2_intv"],
    )
    y_d = _find_value(
        rows,
        lambda r: _same_level(r["main_d"], CONTRAST_D) and r["m1_d"] is None and r["m2_d"] is None and not r["m2_intv"],
    )
    y_d_intv = _find_value(rows, lambda r: _same_level(r["main_d"], CONTRAST_D) and r["m2_intv"])

    if not any(pd.isna(v) for v in [y_dstar, y_d, y_d_intv]):
        # OE = E[Y(d)] - E[Y(d*)]
        # IDE = E[Y(d, M2=intv)] - E[Y(d*)]
        # IIE = E[Y(d)] - E[Y(d, M2=intv)]
        out["mf_intv_OE"] = y_d - y_dstar
        out["mf_intv_IDE"] = y_d_intv - y_dstar
        out["mf_intv_IIE"] = y_d - y_d_intv
    return out

def run_one_rep(
    rep_id,
    output_dir,
    train_kwargs=None,
    result_prefix="medflow",
    model_dir_name="models",
    pse_prefix="pse",
    intv_prefix="intv",
):
    from medflow import train_med, sim_med

    rep_id = int(rep_id)
    rep_dir = os.path.join(output_dir, f"rep_{rep_id}")
    results = {"rep_id": rep_id}
    train_kwargs = dict(train_kwargs or {})

    try:
        # ------------------------------------------------------------------
        # 1) Read CSV data (written by R worker) and prepare for medflow
        # ------------------------------------------------------------------
        data_csv = os.path.join(rep_dir, "data.csv")
        if not os.path.exists(data_csv):
            raise FileNotFoundError(f"{data_csv} not found")

        df = pd.read_csv(data_csv)
        print(f"[Rep {rep_id}] Data loaded: {df.shape[0]} rows, {df.shape[1]} columns")

        # ------------------------------------------------------------------
        # 2) Train cGNF model
        # ------------------------------------------------------------------
        print(f"[Rep {rep_id}] Training medflow model...")
        default_train_kwargs = {
            "path": rep_dir + "/",
            "dataset_name": "data",
            "model_name": model_dir_name,
            "treatment": "D",
            "mediator": ["M1", "M2"],
            "outcome": "Y",
            "confounder": ["C"],
            "cat_var": ["M1", "Y"],
            "nb_epoch": 50000,
            "nb_estop": 50,
            "trn_batch_size": 128,
            "val_batch_size": 2048,
            "learning_rate": 1e-4,
            "seed": rep_id,
            "emb_net": [100, 90, 80, 70, 60],
            "int_net": [60, 50, 40, 30, 20],
        }
        default_train_kwargs.update(train_kwargs)
        train_med(**default_train_kwargs)
        print(f"[Rep {rep_id}] Model training complete.")

        # ------------------------------------------------------------------
        # 3) PSE simulation (intv_med=None)
        # ------------------------------------------------------------------
        print(f"[Rep {rep_id}] Running PSE simulation...")
        pse_result = sim_med(
            path=rep_dir + "/",
            dataset_name="data",
            model_name=model_dir_name,
            cat_list=CAT_LIST,
            intv_med=None,
            n_mce_samples=5000,
            seed=rep_id,
            inv_datafile_name=pse_prefix,
        )
        print(f"[Rep {rep_id}] PSE simulation complete.")

        # ------------------------------------------------------------------
        # 4) Interventional simulation (try both intv_med formats)
        # ------------------------------------------------------------------
        print(f"[Rep {rep_id}] Running interventional simulation...")
        intv_result = sim_med(
            path=rep_dir + "/",
            dataset_name="data",
            model_name=model_dir_name,
            cat_list=CAT_LIST,
            intv_med=["M2=intv"],
            n_mce_samples=5000,
            seed=rep_id + 100000,
            inv_datafile_name=intv_prefix,
        )
        intv_df = _coerce_result_df(intv_result, rep_dir, intv_prefix)
        if intv_df is None or intv_df.empty:
            raise RuntimeError(f"[Rep {rep_id}] Missing interventional medflow summary output.")
        print(f"[Rep {rep_id}] Interventional simulation complete.")

        # ------------------------------------------------------------------
        # 5) Parse results (robust to medflow output schema changes)
        # ------------------------------------------------------------------
        pse_df = _coerce_result_df(pse_result, rep_dir, pse_prefix)
        if pse_df is None or pse_df.empty:
            raise RuntimeError(f"[Rep {rep_id}] Missing PSE medflow summary output.")
        if pse_df is not None:
            print(f"[Rep {rep_id}] PSE columns: {list(pse_df.columns)}")
        if intv_df is not None:
            print(f"[Rep {rep_id}] INTV columns: {list(intv_df.columns)}")

        pse_patterns = {
            "mf_ATE": r"TE\(|\bATE\b",
            "mf_DY": r"PSE.*D\s*[-=~]*>\s*Y|MNDE|direct",
            "mf_DM2Y": r"PSE.*M2\s*[-=~]*>\s*Y|via\s*M2",
            "mf_DM1Y": r"PSE.*M1\s*[-=~]*>\s*Y|via\s*M1",
        }
        intv_patterns = {
            "mf_intv_OE": r"OE\(",
            "mf_intv_IDE": r"IDE\(",
            "mf_intv_IIE": r"IIE\(",
        }

        pse_po = _extract_pse_from_potential_outcomes(pse_df, rep_id)
        intv_po = _extract_intv_from_potential_outcomes(intv_df, rep_id)

        if all(pd.isna(v) for v in pse_po.values()):
            print(f"[Rep {rep_id}] WARNING: PSE potential-outcome parsing failed; retrying with labeled regex extraction.")
            pse_po = _extract_effects(pse_df, pse_patterns, rep_id=rep_id, label="PSE")
        if all(pd.isna(v) for v in intv_po.values()):
            print(f"[Rep {rep_id}] WARNING: INTV potential-outcome parsing failed; retrying with labeled regex extraction.")
            intv_po = _extract_effects(intv_df, intv_patterns, rep_id=rep_id, label="INTV")

        if all(pd.isna(v) for v in pse_po.values()):
            raise RuntimeError(f"[Rep {rep_id}] Unable to extract PSE estimates from medflow output.")
        if all(pd.isna(v) for v in intv_po.values()):
            raise RuntimeError(f"[Rep {rep_id}] Unable to extract interventional estimates from medflow output.")

        results.update(pse_po)
        results.update(intv_po)

    except Exception as e:
        err_dir = os.path.join(output_dir, "error_logs")
        os.makedirs(err_dir, exist_ok=True)
        err_file = os.path.join(err_dir, f"{result_prefix}_rep_{rep_id}_error.log")
        with open(err_file, "w") as f:
            f.write(f"rep_id: {rep_id}\n")
            f.write(f"error: {str(e)}\n\n")
            f.write(traceback.format_exc())
        print(f"[Rep {rep_id}] ERROR: {e}")
        print(f"[Rep {rep_id}] Error log written to {err_file}")

        results["mf_ATE"] = np.nan
        results["mf_DY"] = np.nan
        results["mf_DM2Y"] = np.nan
        results["mf_DM1Y"] = np.nan
        results["mf_intv_OE"] = np.nan
        results["mf_intv_IDE"] = np.nan
        results["mf_intv_IIE"] = np.nan

    # ------------------------------------------------------------------
    # 6) Save results
    # ------------------------------------------------------------------
    rep_result_dir = _rep_result_dir(output_dir)
    results_df = pd.DataFrame([results])
    out_path = os.path.join(rep_result_dir, f"{result_prefix}_rep_{rep_id}.csv")
    results_df.to_csv(out_path, index=False)
    print(f"[Rep {rep_id}] Results saved to {out_path}")

    return results


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python mc_exp4_py_worker.py <task_id> <n_reps> <n_nodes> <n_cores> <output_dir>")
        sys.exit(1)

    task_id    = int(sys.argv[1])
    n_reps     = int(sys.argv[2])
    n_nodes    = int(sys.argv[3])
    n_cores    = int(sys.argv[4])
    output_dir = sys.argv[5]

    # Calculate which replications this task handles
    import math
    from joblib import Parallel, delayed

    reps_per_node = math.ceil(n_reps / n_nodes)
    start_rep = (task_id - 1) * reps_per_node + 1
    end_rep = min(task_id * reps_per_node, n_reps)
    my_reps = list(range(start_rep, end_rep + 1))

    print(f"Task {task_id} handling replications {start_rep} to {end_rep} ({len(my_reps)} total)")
    print(f"Using {n_cores} parallel processes")
    report_runtime_backend("medflow")

    results = Parallel(n_jobs=n_cores)(
        delayed(run_one_rep)(rep_id, output_dir) for rep_id in my_reps
    )
