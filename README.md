# Using Prior Information to Boost Power in Correlation Structure Support Recovery
Focuses on the support recovery problem of correlation matrix structure using the Frequentist Assisted by Bayes approach

# Author Contributions Checklist Form

## Authors
* Ziyang Ding
* David Dunson

## ArXiv Link
https://arxiv.org/abs/2111.11278

## Data

### Abstract
All external datasets come from Cancer DepMap portal, which are open datasets. The Cancer Dependency Map contains an extensive collection of genomics data, including measurements of gene expression, RNAi and CRISPR dependency, and drug sensitivity gathered from over 1000 cancer cell types.

### Availability
All required datasets are incorporated in this github repository. 

### Data Directory
* ./data/coexpressed_genes.csv : This is the correlated gene pairs we use as ground truth for real data experiment as presented in Section 5 of our manuscript
* ./data/depmap_crispr_data.csv: This is this gene expression for breast cancer gene expression dataset used in Section 5, case 1, as the auxiliary dataset
* ./data/depmap_rnai_data.csv: This is this gene expression for breast cancer gene expression dataset used in section 5, case 1, as the testing dataset
* ./data/depmap_breast_gene_expression_data.csv: This is this gene expression for breast cancer gene expression dataset used in Section 5, case 2, as the auxiliary dataset
* ./data/depmap_lung_gene_expression_data.csv: This is this gene expression for lung cancer gene expression dataset used in Section 5, case 2, as the testing dataset

## Code
The root folder contains all necessary files to implement all the algorithms, figure and table generations, and simulation data generations. There are 3 major kinds of codde
* .R files: The R scirpts include algorithms, figure generators, and simulated data generators
* .Rmd files: The Rmd scipts doesn’t carry any such functionality but calls functions in R scripts to form reproducibility files. All results in the the paper are generated using Rmd files.
* .RCpp file: There is a single RCpp file intended to speed up computation.

### Code File Explanation
Below is a summary of each code file:
* [Dependency.R] (Dependency.R): This file contains all the packages and dependency specification. It will be automatically run everytime the reader runs any Rmd file
* [Utils.R] (Utils.R): This file contains all algorithms, figure and table generations, and simulation data generations.
* [fun_sim_idpt.R] (fun_sim_idpt.R): This file runs simulation for the independent FAB correlation structure test as in Secion 4.1 in the paper.
* [fun_sim_boot.R] (fun_sim_boot.R): This file runs simulation for the bootstrap FAB correlation structure test as in Secion 4.2 in the paper.
* [software.R] (software.R): This file contains the simulation function used in the real_data_experiment.rmd. 
* [IndividualExperiments.Rmd] (IndividualExperiments.Rmd): This Rmd file reproduces all individual simulation experiments in Section 4 in the paper.
* [MassiveExperiments.Rmd] (MassiveExperiments.Rmd): This Rmd file reproduces all composite simulation experiments. The composite simulation experiments means simulation experiments that requires generatation of the dataset for multiple times and average their final results.
* [real_data_experiments] (real_data_experiments): This Rmd file reproduces all real data experiments in Section 5 in the paper.
* [Cov_struct_func.cpp] (Cov_struct_func.cpp): an RCpp file that speeds up computation for the bootstrap FAB.

### Dependency and Version Control
We uses R version 4.0.5. The dependency packages version is listed in [version_control.txt] (version_control.txt)

## Reproducing Workflow
Our code reproduces all numbers, tables, and figures in the paper. The following list summarizes workflows to reproduce the required results.
* To reproduce individual simulation section results: open individualExperiments.Rmd, simply run all chunks.
* To reproduce real data section results: open real_data_experiments.Rmd, simply run all chunks
* To reproduce table and bootstrap grid figure in our paper: open MassiveExperiments.Rmd, simply run all chunks

**Please note**: please restart R session and clear all output everytime run another Rmd file. Otherwise, results could be slightly different.
