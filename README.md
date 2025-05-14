# Causal Mediation Analysis with Multiple Mediators: A Simulation Approach

This repository contains the replication files for the paper **"Causal Mediation Analysis with Multiple Mediators: A Simulation Approach"**.

## Replication Data

The paper reanalyzes the 2003 U.S. birth certificate data used in VanderWeele et al. (2014) and the media framing experiment data from Brader et al. (2008), building on related work by Imai et al. (2013), Zhou and Yamamoto (2019), and Wodtke and Zhou (forthcoming).

- The **2003 U.S. birth certificate data** can be downloaded here: [Download Data](https://data.nber.org/lbid/2003/linkco2003us_den.csv.zip)  
- The **media framing experiment data** can be downloaded here: [Download Data](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/OU6D17)


## Files and descriptions:

- **`Table_1.R`**: Performs parametric causal mediation analysis to replicate the results in Table 1.

### Folder: `Table_2`
Contains scripts to replicate results in Table 2.

- **`clean.R`**: Cleans the raw dataset for analysis.
- **`UMNNs.py`**: Performs causal mediation analysis using normalizing flows and UMNNs to replicate the results in the first column.
- **`Parametric.R`**: Performs parametric causal mediation analysis to replicate the results in the second column.


### Folder: `Figure_4`
Contains scripts for replicating the plots in Figure 4.

- **`Dequantization.py`**: Generates the dequantized data from the original data with normalizing flows .
- **`Normalizing.py`**: Transforms the original data into a standard Gaussian distribution with normalizing flows.
- **`Plot.py`**: Plots the final figure.

## Instructions:

1. Install the [`MedSim`](https://github.com/causalMedAnalysis/causalMedR) package in R for parametric causal mediation analysis.
2. Install the [`MedFlow`](https://github.com/JesseZhou-1/medflow) package in Python for causal mediation analysis using normalizing flows and UMNNs.
