# Causal Mediation Analysis with Multiple Mediators: A Simulation Approach

This repository contains the replication files for the paper **"Causal Mediation Analysis with Multiple Mediators: A Simulation Approach"**.

The repository includes code and data for:

1. the empirical application on media framing and immigration attitudes in the `Immigration` folder;
2. the empirical application on prenatal care and preterm birth in the `Preterm Birth` folder; and
3. the five Monte Carlo simulation experiments in the `Monte Carlo` folder.

## Replication Data

The paper reanalyzes the 2003 U.S. birth certificate data used in VanderWeele et al. (2014) and the media framing experiment data from Brader et al. (2008), building on related work by Imai et al. (2013), Zhou and Yamamoto (2019), and Wodtke and Zhou (forthcoming).

- The **2003 U.S. birth certificate data** can be downloaded here: [Download Data](https://data.nber.org/lbid/2003/linkco2003us_den.csv.zip)  
- The **media framing experiment data** can be downloaded here: [Download Data](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/OU6D17)

## Code and descriptions:

### Folder: `Immigration`

- **`Table_6.R`**: replicates the empirical results for the media-framing application reported in Table 6 of the manuscript.`

### Folder: `Preterm Birth`

#### Folder: `Table 7`

- **`Clean.R`**: cleans and prepares the linked birth-certificate data.
- **`UMNNs.py`**: runs the normalizing-flow / UMNN specification.
- **`Parametric.R`**: runs the parametric simulation estimator.

#### Folder: `Figure_4`

- **`Dequantization.py`**: generates the dequantized variables used in the diagnostic plots.
- **`Normalizing.py`**: applies the trained normalizing flows.
- **`Plot.py`**: produces the final diagnostic figure.

### Folder: `Monte Carlo`

#### Floder: `auxiliary`

This subfolder contains shared estimator implementations and support functions used across experiments, including:

* `medsim.R`: simulation-based estimator;
* `linpath.R`: linear-model estimator for path-specific effects;
* `ipwpath.R` and `ipwmed.R`: IPW estimators for path-specific and related effects;
* `ipwvent.R`: IPW estimator for interventional effects;
* `pathimp.R`: regression-imputation estimator;
* `rwrlite.R` and `rwrmed.R`: regression-with-residuals code;
* `dmlmed.R`, `dmlpath.R`, and `dmlpath_custom.R`: DML-based estimators used in the nonlinear experiments.

#### Experiment folders

`exp1`-`exp5` contain scripts used to replicate Experiment 1-5. Within these folders:

* `run_exp*.R` scripts are stand-alone R drivers for local checking, development, and smaller-scale runs;
* `mc_*_worker.*` scripts are worker scripts used for larger Monte Carlo runs;
* `aggregate_results_*.R` scripts combine replication-level outputs and compute bias and RMSE summaries;
* `run_mc_*.sbatch and run_mc_*.sh` scripts are SLURM/HPC launch scripts used for cluster execution.
  
Experiment 4 combines R-based estimators with a Python MedFlow workflow. Experiment 5 compares MedFlow hyperparameter variants using the Experiment 4 design.

### Attribution and Provenance
Parts of the auxiliary function code build on earlier causal mediation replication code associated with:
* the repFiles repository for Causal Mediation Analysis: https://github.com/causalMedAnalysis/repFiles/tree/50e575f284b2312d7c4189f33a69472a53e1ca1c
* the accompanying book: https://www.cambridge.org/us/universitypress/subjects/social-science-research-methods/quantitative-methods/causal-mediation-analysis
  
The regression-imputation code also relies on a modified version of the paths package. Relevant upstream sources are:
* modified fork used for this project: https://github.com/JesseZhou-1/paths
* original upstream repository: https://github.com/xiangzhou09/paths

The normalizing-flow code requires the installation of the [`MedFlow`](https://github.com/JesseZhou-1/medflow) package in Python.
