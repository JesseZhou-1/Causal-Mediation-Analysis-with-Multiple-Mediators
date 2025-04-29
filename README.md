# Causal Mediation Analysis with Multiple Mediators: A Simulation Approach

This repository contains the replication files for the paper **"Causal Mediation Analysis with Multiple Mediators: A Simulation Approach"**.

## Replication Data

The paper reanalyzes the 2003 U.S. birth certificate data used in VanderWeele et al. (2014) and the media framing experiment data from Brader et al. (2008), building on related work by Imai et al. (2013), Zhou and Yamamoto (2019), and Wodtke and Zhou (forthcoming).

- The **2003 U.S. birth certificate data** can be downloaded here: [Download Data](https://data.nber.org/lbid/2003/linkco2003us_den.csv.zip)  
- The **media framing experiment data** can be downloaded here: [Download Data](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/OU6D17)


## Files and Descriptions:

- **`clean.R`**: Script to clean the raw dataset for analysis.
- **`deq.py`**: Generates the dequantized data using cGNF for Figure 6.
- **`normalizing.py`**: Transforms the data into a standard Gaussian distribution using cGNF for Figure 6.
- **`transformation_plot.py`**: Combines the plots of the original, dequantized, and transformed data for Figure 6.
- **`medsim_replication.R`**: Performs parametric causal mediation analysis to replicate results in Table 1.
- **`medsimGNF_replication.R`**: Executes causal mediation analysis using cGNF to replicate results in Table 2.
- **`sensitivity.py`**: Conducts sensitivity analysis using cGNF to replicate results in Table 2.

## Instructions:

1. Install the [`MedSim`](https://github.com/JesseZhou-1/medsim) package in R for parametric causal mediation analysis.
2. Install the [`MedFlow`](https://github.com/JesseZhou-1/medflow) package in Python for causal mediation analysis using cGNF.
