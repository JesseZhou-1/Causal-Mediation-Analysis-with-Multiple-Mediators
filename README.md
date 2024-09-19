# Causal Mediation Analysis with Multiple Mediators: A Simulation Approach: Replication Files

This repository contains the replication files for the paper **"Causal Mediation Analysis with Multiple Mediators: A Simulation Approach"**.

## Files and Descriptions:

- **`clean.R`**: Script to clean the raw dataset for analysis.
- **`deq.py`**: Generates the dequantized data using cGNF for Figure 6.
- **`normalizing.py`**: Transforms the data into a standard Gaussian distribution using cGNF for Figure 6.
- **`transformation_plot.py`**: Combines the plots of the original, dequantized, and transformed data for Figure 6.
- **`medsim_replication.R`**: Performs parametric causal mediation analysis to replicate results in Table 1.
- **`medsimGNF_replication.R`**: Executes causal mediation analysis using cGNF to replicate results in Table 2.
- **`sensitivity.py`**: Conducts sensitivity analysis using cGNF to replicate results in Table 2.

## Instructions:

1. Install the `medsim` package in R for parametric causal mediation analysis.
2. Install the `medsimGNF` package in Python for causal mediation analysis using cGNF.
