# Uncovering the link between gut microbiome, highly-processed food consumption and diet quality through bioinformatics methods

This repository contains all scripts and auxiliary data used to analyse metagenomic data from two cohorts, perform differential abundance analyses (DAA) and build all classifiers, as well as regarding the validation of one of these ML models and the figures shown in this MSc Thesis memory.

## Table of contents
* [Data directory](#data-directory)
* [Files for heatmap construction](#files-for-heatmap-construction)
* [R scripts](#r-scripts)
* [Taxonomic assignment](#taxonomic-assignment)

## Data directory

This directory contains all files used for the training of all Machine Learning classifiers, for the generation of both heatmaps shown in the TFM Thesis and for the plotting of all figures shown in the same document.

- `data_clr_txt`: text file

## Files for heatmap construction

The `heatmap_files` directory contains all data necessary for the construction of heatmaps shown in the MSc Thesis (Figures 4 and 7):
- One directory called `additional_files_hpf_classification`. Three files are stored here:
    - `clr_abundances_for_heatmap.csv`: this file contains CLR-transformed relative abundances for the 20 most important taxa according to the most precise ML classifier based on HPF consumption, where each row is a sample and each column is a genus. Samples are ordered so that those belonging to the "Low HPF consumption" group appear first, followed by those from the "High HPF consumption" group.
    - `groups_for_heatmap_genuses.csv`: this file consists in 20 rows and 2 columns: rows are for each one of the 20 most important taxa according to the best classifier based on this classification of the study cohort, while the first column are these taxa and the second one indicated the HPF consumption group in which they are enriched.
    - `groups_for_heatmap.csv`: this file indicates which HPF consumption group each sample of the study cohort belongs to. There are therefore as many rows as samples there are in this cohort, while the second column details the group each individual pertains to.
- One directory called `additional_files_hei_classification`. Three files are stored here:
    - `clr_abundances_for_heatmap_hei.csv`: this file contains CLR-transformed relative abundances for the 20 most important taxa according to the most precise ML classifier based on HEI, where each row is a sample and each column is a genus. Samples are ordered so that those belonging to the "Good HEI" group appear first, followed by those from the "Poor HEI" group.
    - `groups_for_heatmap_genuses_hei.csv`: this file consists in 20 rows and 2 columns: rows are for each one of the 20 most important taxa according to the best classifier based on this classification of the study cohort, while the first column are these taxa and the second one indicated the HEI group in which they are enriched.
    - `groups_for_heatmap_hei.csv`: this file indicates which HEI group each sample of the study cohort belongs to. There are therefore as many rows as samples there are in this cohort, while the second column details the group each individual pertains to.

## R scripts

The `R_scripts` directory contains:
- `alpha_and_beta_diversity.R`:
- `best_ML_model_study_cohort_HEI.R`:
- `boxplots.R`:
- `daa.R`:
- `importance_barplots.R`:
- `metaphlan_to_R.R`:
- `ML_models_HEI.R`:
- `ML_models_HPF.R`:
- `population_tables.R`:
- `qiime_to_R.R`:
- `quantile_graphs.R`:
- `roc_curves.R`:
- `validation.R`:

## Taxonomic assignment

The `taxonomic_assignment` directory contains:
- `16s_data.sh`: bash script used for the processing of raw reads using QIIME2 and their taxonomic assignment, until obtaining all necessary files for the generation of phyloseq objects that will be needed for ML classifier construction based on HPF consumption and HEI.
- `wgs_data.sh`: bash script used for the processing of raw reads using metaWRAP and MetaPhlAn4 and their taxonomic assignment, until obtaining the relative abundances table necessary for the generation of phyloseq objects that will be needed for ML classifier construction based on HEI.

