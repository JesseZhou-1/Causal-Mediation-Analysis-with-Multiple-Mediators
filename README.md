# Causal Mediation Analysis with Multiple Mediators: A Simulation Approach

This repository contains replication files for the paper **"Causal Mediation
Analysis with Multiple Mediators: A Simulation Approach."** It includes:

1. the media-framing and immigration-attitudes application in `Immigration`;
2. the prenatal-care and preterm-birth application in `Preterm Birth`; and
3. five Monte Carlo experiments in `Monte Carlo`.

## Replication Data

The paper reanalyzes the 2003 U.S. birth-certificate data used by VanderWeele
et al. (2014) and the media-framing experiment from Brader et al. (2008).

- [2003 U.S. birth-certificate data](https://data.nber.org/lbid/2003/linkco2003us_den.csv.zip)
- [Media-framing experiment data](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/OU6D17)

Data files are not included in this repository. Download each file and place it
beside the corresponding cleaning or analysis script before running the code.

## Software

The R scripts use standard CRAN packages, including `MASS`, `SuperLearner`,
`caret`, `dplyr`, `nnet`, `purrr`, `ranger`, `readr`, `stringr`,
and `tidyr`.

The regression-imputation code requires the modified `paths` package used for
this project:

```r
devtools::install_github("JesseZhou-1/paths")
```

The normalizing-flow analyses require
[`MedFlow`](https://github.com/JesseZhou-1/medflow):

```bash
python -m pip install --upgrade git+https://github.com/JesseZhou-1/medflow.git
```

The RQS scripts additionally require the spline-enabled
[`cGNF-spline`](https://github.com/JesseZhou-1/cGNF-spline) backend. It
provides both the UMNN and RQS normalizers, so it can be used for all MedFlow
scripts in this repository. Replace the original `cGNF` backend as follows:

```bash
python -m pip uninstall -y cGNF cGNF-spline
python -m pip install git+https://github.com/JesseZhou-1/cGNF-spline.git
```

The Slurm scripts contain resource requests and module commands for the
University of Chicago Midway cluster. Users on other systems should adjust
those directives to their computing environment.

## Repository Structure

### `Immigration`

- `Table_6.R` reproduces the media-framing results reported in Table 6.

### `Preterm Birth/Table_7`

- `Clean.R` cleans and prepares the linked birth-certificate data.
- `UMNNs.py` fits the MedFlow-UMNN specification.
- `RQS.py` fits the MedFlow-RQS specification.
- `Parametric.R` fits the parametric simulation estimator.

`RQS.py` defaults to the reported RQS architecture: embedding layers of
`[100, 90, 80, 70, 60]`, spline-parameter layers of
`[60, 50, 40, 30, 20]`, 256 bins, and a bound of 30 on the standardized
scale. It uses seed `250470397` and 100,000 Monte Carlo draws. These settings
can be changed through `MEDFLOW_SEED`, `MEDFLOW_N_MCE`, and
`MEDFLOW_DATA_DIR`. Training or either estimation stage can be skipped with
`MEDFLOW_RUN_TRAIN=0`, `MEDFLOW_RUN_PATH=0`, or `MEDFLOW_RUN_INTV=0`.

After downloading `linkco2003us_den.csv`, run:

```bash
cd "Preterm Birth/Table_7"
Rscript Clean.R
python RQS.py
```

### `Preterm Birth/Figure_4`

- `Dequantization.py` generates the dequantized variables.
- `Normalizing.py` applies the trained UMNN normalizing flows.
- `Plot.py` produces the diagnostic figure.

As stated in the manuscript, Figure 4 illustrates the UMNN model only.

### `Monte Carlo/auxiliary`

This folder contains shared estimator implementations and support functions:

- `medsim.R`: simulation-based estimator;
- `linpath.R`: linear-model estimator for path-specific effects;
- `ipwpath.R`, `ipwmed.R`, and `ipwvent.R`: IPW estimators;
- `pathimp.R`: regression-imputation estimator;
- `rwrlite.R` and `rwrmed.R`: regression-with-residuals code; and
- `dmlmed.R`, `dmlpath.R`, and `dmlpath_custom.R`: DML estimators.

### `Monte Carlo/exp1` through `Monte Carlo/exp3`

Each folder contains a stand-alone R driver, an HPC worker and launcher, and an
aggregation script. The local `run_exp*.R` files are useful for checking the
workflow on a smaller scale; the `run_mc_*.sbatch` files reproduce the full
Monte Carlo runs.

### `Monte Carlo/exp4`

Experiment 4 first runs the R worker, which generates each replication's data
and computes the R-based estimators. Separate Python workers then fit the
MedFlow-UMNN and MedFlow-RQS models. The RQS worker uses the settings reported
in the manuscript: a training batch size of 128, validation batch size of
2,048, learning rate of 0.0001, 256 bins, bound 7, and spline-parameter layers
of `[60, 50, 40, 30, 20]`.

On Slurm, run:

```bash
cd "Monte Carlo/exp4"
sbatch run_mc_exp4_r.sbatch

# Set this to the directory created by the R job.
EXP4_OUTPUT_DIR=/path/to/results_exp4_JOBID

sbatch run_mc_exp4_py.sbatch "$EXP4_OUTPUT_DIR"
sbatch run_mc_exp4_py_spline.sbatch "$EXP4_OUTPUT_DIR"
```

The two Python workers may use the same Experiment 4 directory. UMNN outputs
use `models` and `medflow_rep_*.csv`; RQS outputs use `models_spline` and
`medflow_spline_rep_*.csv`, so neither run overwrites the other.

### `Monte Carlo/exp5`

Experiment 5 reuses the first 240 data sets from Experiment 4 and varies model
capacity, learning rate, and training batch size. The original launchers run
the UMNN variants. The `*_spline_*` launchers run the corresponding RQS
variants. For RQS, the architecture variants jointly change the embedding
network, spline-parameter network, and number of bins.

```bash
cd "Monte Carlo/exp5"

EXP4_OUTPUT_DIR=/path/to/results_exp4_JOBID
EXP5_RQS_DIR=/path/to/results_exp5/spline

sbatch run_mc_exp5_spline_arch.sbatch  "$EXP4_OUTPUT_DIR" "$EXP5_RQS_DIR" up
sbatch run_mc_exp5_spline_arch.sbatch  "$EXP4_OUTPUT_DIR" "$EXP5_RQS_DIR" down
sbatch run_mc_exp5_spline_lr.sbatch    "$EXP4_OUTPUT_DIR" "$EXP5_RQS_DIR" up
sbatch run_mc_exp5_spline_lr.sbatch    "$EXP4_OUTPUT_DIR" "$EXP5_RQS_DIR" down
sbatch run_mc_exp5_spline_batch.sbatch "$EXP4_OUTPUT_DIR" "$EXP5_RQS_DIR" up
sbatch run_mc_exp5_spline_batch.sbatch "$EXP4_OUTPUT_DIR" "$EXP5_RQS_DIR" down
```

Each launcher creates its own variant subdirectory and aggregates that
variant. After all six jobs finish, combine the RQS summaries with:

```bash
Rscript compare_results_exp5.R \
  "$EXP4_OUTPUT_DIR" "$EXP5_RQS_DIR" 240 medflow_spline
```

The optional fourth argument to `aggregate_results_exp5_variant.R` and
`compare_results_exp5.R` selects the result prefix. It defaults to
`medflow` for UMNN and should be `medflow_spline` for RQS.

## Attribution and Provenance

Parts of the auxiliary code build on the replication files for the book
[*Causal Mediation Analysis*](https://www.cambridge.org/us/universitypress/subjects/social-science-research-methods/quantitative-methods/causal-mediation-analysis),
available in the
[`causalMedAnalysis/repFiles` repository](https://github.com/causalMedAnalysis/repFiles/tree/50e575f284b2312d7c4189f33a69472a53e1ca1c).

The regression-imputation implementation uses a modified fork of the
[`paths` package](https://github.com/JesseZhou-1/paths). The
[original package repository](https://github.com/xiangzhou09/paths) is
maintained by Xiang Zhou.

The normalizing-flow implementation uses
[`MedFlow`](https://github.com/JesseZhou-1/medflow),
[`cGNF`](https://github.com/cGNF-Dev/cGNF), and the
[`cGNF-spline`](https://github.com/JesseZhou-1/cGNF-spline) extension.
