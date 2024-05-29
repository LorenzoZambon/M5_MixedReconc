# M5 - Mixed reconciliation

Code for reproducing the results of the experiments on the M5 dataset of the paper *Probabilistic reconciliation of mixed-type hierarchical time series* [@zambon2024mixed], 
published at UAI 2024 (the 40th Conference on Uncertainty in Artificial Intelligence).
We refer to the paper and to the vignette *Reconciliation of M5 hierarchy with mixed-type forecasts* of the R package `bayesRecon` for all the details on the experiments.

## Requirements

This is the list of the required packages, available on CRAN:
* bayesRecon
* m5
* parallel
* doSNOW
* foreach
* smooth
* tictoc

## Usage

Download the repository and run the file main.R.
The folders where data and results are saved can be set in main.R.

## Cite Us
```
@inproceedings{
zambon2024mixed,
title={Probabilistic reconciliation of mixed-type hierarchical time series},
author={Lorenzo Zambon and Dario Azzimonti and Nicolò Rubattu and Giorgio Corani},
booktitle={The 40th Conference on Uncertainty in Artificial Intelligence},
year={2024}
}
```
