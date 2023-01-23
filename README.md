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
- `alpha_and_beta_diversity.R`: R script describing all procedures needed to plot calculate alpha diversity measures and plot several indices for both cohorts and classifications.
- `best_ML_model_study_cohort_HPF.R`: R script describing how the most precise model for the HPF-based classification was built. Feature importances are calculated as well. This script also includes the code needed to draw the heatmap shown in Figure 4 of this MSc Thesis, relative to the HPF classificafion of the study cohort.
- `best_ML_model_study_cohort_HEI.R`: R script describing how the most precise model for the HEI-based classification was built. Feature importances are calculated as well. This script also includes the code needed to draw the heatmap shown in Figure 7 of this MSc Thesis, relative to the HEI classificafion of the study cohort.
- `boxplots.R`: R script indicating in detail how boxplots showing taxa being significantly different in their abundances between HPF- or HEI-based groups were made. It is therefore indicated how those individuals with type 2 diabetes family history were highlighted, as well as how the Wilcoxon rank sum test was used to find these differentially abundant taxa between groups.
- `daa.R`: R script describing how differential abundance analyses were performed for both classifications of the study cohort and the validation one grouped based on their HEI. 
- `importance_barplots.R`: R script describing how supplementary figures, consisting in barplots and showing the most important taxa according to each one of our ML classifiers, were made.
- `metaphlan_to_R.R`: R script indicating how, once the relative abundances table from the validation cohort was obtained after executing MetaPhlAn4, these abundaces were imported to R for the creation of all phyloseq objects necessary for diversity analyses, DAA and ML classifier validation.
- `ML_models_HEI.R`: R script indicating in detail how all 6 ML models for the HEI classification of the study cohort were built.
- `ML_models_HPF.R`: R script indicating in detail how all 6 ML models for the HPF classification of the study cohort were built.
- `population_tables.R`: R script describing how Tables 1-3, which showed some demographic and additional data regarding both cohorts, were created using the gtsummary R package.
- `qiime_to_R.R`: R script indicating how, once all QIIME2-related procedures for the study cohort were finished, the resulting artifacts were used for the creation of phyloseq objects in R necessary for diversity analyses, DAA and ML classifier construction.
- `quantile_graphs.R`: R script describing how supplementary figures, showing boxplots of some differently abundant taxa vs quantiles made using the validation cohort based on their levels of several glycaemic traits (HOMA-IR, HbA1c (%) and adiponectin) and HEI, were obtained. The identification of these taxa differentially abundant among these quantiles via the Kruskal-Wallis rank sum test is also detailed in this script.
- `roc_curves.R`: R script detailing how ROC curves for all ML classifiers regarding both classifications of the study cohort were made.
- `validation.R`: R script describing how validation of the best ML classifier trained on the study cohort classified based on their HEI was performed using the relative abundances table from the validation cohort, along with the plotting of the ROC curve linked to this validation.

## Taxonomic assignment

The `taxonomic_assignment` directory contains:
- `16s_data.sh`: bash script used for the processing of raw reads using QIIME2 and their taxonomic assignment, until obtaining all necessary files for the generation of phyloseq objects that will be needed for ML classifier construction based on HPF consumption and HEI.
- `wgs_data.sh`: bash script used for the processing of raw reads using metaWRAP and MetaPhlAn4 and their taxonomic assignment, until obtaining the relative abundances table necessary for the generation of phyloseq objects that will be needed for ML classifier construction based on HEI.

