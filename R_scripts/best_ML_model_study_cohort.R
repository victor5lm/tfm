# Scripts for the generation of the most accurate ML classifier for the HPF classification of the Study cohort:

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(phyloseq)
library(viridis)
library(RColorBrewer)
library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(purrr)
library(mikropml)

# Classifier construction

library(caret)
library(pROC)
library(ggpubr)
library(tidyverse)
library(viridis)

if( !require("mikropml") ){
    install.packages("mikropml")
    library(mikropml)
}

data_complete <- read.table('data_clr_complete.tsv',
                      header = TRUE, sep = '\t',# skip = 1,
                      row.names = 1)
ntrees<-seq(1,100)
tuning_rf<-list(mtry=seq(1,70))
originals <- list()

for (ntree in ntrees){print(ntree)
    set.seed(2019)
    results <- run_ml(data_complete, "rf",
    outcome_colname = 'Grupo_HPF',
    cross_val=caret::trainControl(method = "LOOCV"),
    training_frac = 0.80, seed=2019,
    calculate_performance = TRUE, ntree=ntree,
    hyperparameters = tuning_rf,
    find_feature_importance = TRUE)
                      
    key<-toString(ntree)
    originals[[key]]<-results
}

# Heatmap for the 20 most important microbial taxa identified by the best ML model trained after the HPF classification of the study cohort

# Heatmap for the 20 most important microbial taxa identified by the best ML model trained after the HEI classification of the study cohort
