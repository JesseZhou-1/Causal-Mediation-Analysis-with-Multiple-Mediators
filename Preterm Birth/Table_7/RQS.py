#!/usr/bin/env python3
"""Fit the MedFlow-RQS model used in the preterm-birth application.

Run ``Clean.R`` first so that ``cleaned_natl2003_bin.csv`` is available.
By default, data and outputs are read from and written to this script's
directory. Set ``MEDFLOW_DATA_DIR`` to use another directory.
"""

from __future__ import annotations

import collections.abc
import inspect
import os
import sys
import time
from pathlib import Path

import numpy as np

collections.Iterable = collections.abc.Iterable
np.set_printoptions(precision=3, suppress=None)
os.environ["CUDA_VISIBLE_DEVICES"] = os.environ.get("CUDA_VISIBLE_DEVICES", "0")

from medflow import sim_med, train_med


def env_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def report_backend() -> None:
    import cGNF
    from cGNF import train as cgnf_train
    from cGNF.GNF_Modules.Normalizers import SplineNormalizer

    print("medflow path:", Path(sys.modules["medflow"].__file__).resolve())
    print("cGNF path:", Path(cGNF.__file__).resolve())
    print("cGNF.train signature:", inspect.signature(cgnf_train))
    print("SplineNormalizer:", SplineNormalizer.__name__)


def main() -> None:
    start_time = time.time()
    report_backend()

    data_dir = Path(
        os.environ.get("MEDFLOW_DATA_DIR", str(Path(__file__).resolve().parent))
    ).resolve()
    dataset_name = os.environ.get("MEDFLOW_DATASET", "cleaned_natl2003_bin")
    seed = int(os.environ.get("MEDFLOW_SEED", "250470397"))
    n_mce_samples = int(os.environ.get("MEDFLOW_N_MCE", "100000"))
    run_train = env_bool("MEDFLOW_RUN_TRAIN", True)
    run_path = env_bool("MEDFLOW_RUN_PATH", True)
    run_intv = env_bool("MEDFLOW_RUN_INTV", True)

    data_dir.mkdir(parents=True, exist_ok=True)
    dataset_csv = data_dir / f"{dataset_name}.csv"
    if not dataset_csv.exists():
        raise FileNotFoundError(f"Dataset CSV not found: {dataset_csv}")

    model_name = str(data_dir / f"seed_{seed}_spline")
    path = str(data_dir) + os.sep

    print("data_dir:", data_dir)
    print("dataset_name:", dataset_name)
    print("model_name:", model_name)
    print("seed:", seed)
    print("n_mce_samples:", n_mce_samples)

    if run_train:
        train_med(
            path=path,
            dataset_name=dataset_name,
            treatment="a",
            confounder=["c1", "c2", "c3", "c4"],
            mediator=["l", "m"],
            outcome="y",
            test_size=0.2,
            cat_var=["a", "l", "m", "y", "c1", "c2", "c3", "c4"],
            sens_corr=None,
            seed_split=seed,
            model_name=model_name,
            trn_batch_size=128,
            val_batch_size=2048,
            learning_rate=1e-4,
            seed=seed,
            nb_epoch=50000,
            nb_estop=500,
            val_freq=1,
            emb_net=[100, 90, 80, 70, 60],
            int_net=[60, 50, 40, 30, 20],
            norm_type="spline",
            num_bins=256,
            bound=30,
            spline_hidden_dims=(60, 50, 40, 30, 20),
        )

    if not Path(model_name).exists():
        raise FileNotFoundError(f"Model directory not found: {model_name}")

    if run_path:
        sim_med(
            path=path,
            dataset_name=dataset_name,
            model_name=model_name,
            n_mce_samples=n_mce_samples,
            seed=seed,
            inv_datafile_name=f"{seed}_spline_path_{n_mce_samples}",
            cat_list=[0, 1],
            moderator=None,
        )

    if run_intv:
        sim_med(
            path=path,
            dataset_name=dataset_name,
            model_name=model_name,
            n_mce_samples=n_mce_samples,
            intv_med=["m=intv"],
            seed=seed,
            inv_datafile_name=f"{seed}_spline_intv_{n_mce_samples}",
            cat_list=[0, 1],
            moderator=None,
        )

    print("elapsed_seconds:", time.time() - start_time)


if __name__ == "__main__":
    main()
