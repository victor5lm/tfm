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

#Classifier construction

#Variable importances

#Heatmap with the 20 most important features according to this ML classifier and based on the HEI classification of the Study cohort

