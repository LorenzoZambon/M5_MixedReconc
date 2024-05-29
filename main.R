# Set path to current directory
dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)
# Clean environment
rm(list = ls()); gc(); 

library(bayesRecon)
library(m5)
library(parallel)
library(doSNOW)
library(foreach)
library(smooth)
library(tictoc)

source("utils.R")
source("base_forecasts.R")
source("reconcile_forecasts.R")
source("compute_scores.R")
source("boxplot.R")          

set.seed(42)

##################################
# Set parameters

# maximum number of CPUs to use for parallel computation of the base forecasts:
n_cpu = detectCores() - 2  

# methods to compare in the boxplots: 
# choose between ("gauss", "gauss_trunc", "mixed_cond", "TD_cond")  
BP_met = c("gauss", "mixed_cond", "TD_cond")  

data_path    = "../data/"     # where data are saved 
results_path = "../results/"  # where results are saved 

##################################
# Parameters for computing the base forecasts, the reconciled forecasts and the scores
# Do not change if you want to reproduce the results of the paper

STORES = c("CA_1", "CA_2", "CA_3", "CA_4", 
           "WI_1", "WI_2", "WI_3",
           "TX_1", "TX_2", "TX_3")
CAT = c("HOBBIES", "HOUSEHOLD", "FOODS")
DEPT = c("HOBBIES_1", "HOBBIES_2", "HOUSEHOLD_1", "HOUSEHOLD_2", 
         "FOODS_1", "FOODS_2", "FOODS_3")

h_list = 1:14     
n_samples_b = 1e5     # number of samples for the base forecasts of the bottom
N_samples_IS = 1e5    # number of samples for mixedCond
N_samples_TD = 1e4    # number of samples for TDcond
alpha_mis = 0.1       # for Mean Interval Score (e.g. alpha = 0.1 --> 90% intervals)
n_samp_gtrunc = 1e4   # for computing scores of the truncated Gaussian

##################################
# If data and results folder are not present, create them

if (!dir.exists(data_path)) {
  dir.create(data_path)
}
if (!dir.exists(results_path)) {
  dir.create(results_path)
}

##################################
# Download data if not already in data_path

if (!file.exists(paste0(data_path, "sales_test_validation.csv"))) {
  m5::m5_download(data_path)
}

STORES = STORES[1:6]

##################################
# Compute base forecasts

for (STORE in STORES) {
  print(paste0("Computing base forecasts of store ", STORE))
  base_fc_store(STORE, h_list, n_samples_b, 
                data_path, results_path, n_cpu)
}

##################################
# Reconcile base forecasts

for (STORE in STORES) {
  rec_fc_store(STORE, h_list, results_path, 
               N_samples_IS, N_samples_TD)
}

##################################
# Compute scores

for (STORE in STORES) {
  print(paste0("Computing scores for store ", STORE))
  compute_scores(STORE, h_list, alpha_mis, n_samp_gtrunc, results_path)
}

##################################
# Compute the skill scores and produces the boxplots
# Also computes and save the mean skill scores
produce_boxplots(STORES, h_list, BP_met, results_path)

##################################
# Print the mean skill scores
mean_SS = readRDS(paste0(results_path, "mean_skill_scores.rds"))
print("Skill scores for the upper time series: ")
print(mean_SS$upper)
print("Skill scores for the bottom time series: ")
print(mean_SS$bottom)












